#include "mt_data.h"
#include <cassert>
#include <map>
#include "ncl/ncl.h"
#include "ncl/nxsblock.h"
#include "ncl/nxspublicblocks.h"
#include "ncl/nxscxxdiscretematrix.h"
#include "ncl/nxsmultiformat.h"
#include "INIReader.h"
#include "mt_instance.h"
#include "mt_ini_options.h"
#include "mt_char_model.h"
#include "mt_tree.h"
#include "search.h"
#include "pattern_class.h"
using namespace std;

extern bool gQuietMode;
extern long gStrictLevel;
extern bool gValidateInternals;

namespace mt {
    class NCL2MT {
        public:
            void processTree(std::ostream * os,
                unsigned numTaxa,
                const NxsCharactersBlock * charsBlock,
                const NxsDiscreteDatatypeMapper * dataMapper,
                const NxsCDiscreteStateSet ** compressedMatrix,
                const vector<double> & patternWeights,
                const vector<int> &origToCompressed,
                const std::map<unsigned, std::set<unsigned> > & numStates2PatternIndexSet,
                const NxsSimpleTree & tree,
                const ModelDescription & md,
                const INIBasedSettings & ibs);
            INIBasedSettings configureBasedOnINI(MTInstance &,
                                                   ::INIReader & iniReader,
                                                   std::ostream & err);
    };
}

bool preDataINICheck(INIReader & iniReader, std::ostream & err);
int processFilepath(
    const char * filename, // name of the file to be read
    std::ostream * os, // output stream to use (NULL for no output). Not that cerr is used to report errors.
    MultiFormatReader::DataFormatType fmt, // enum indicating the file format to expect.
    INIReader & iniReader);
int readFilepathAsNEXUS(const char *filename, MultiFormatReader::DataFormatType fmt, INIReader & iniReader);
int readFilepathAsNEXUS(const char *filename, MultiFormatReader::DataFormatType fmt, INIReader & iniReader);
int readFilesListedIsFile(const char *masterFilepath, MultiFormatReader::DataFormatType fmt, INIReader & iniReader);
void printHelp(std::ostream & out);
int do_main(int argc, char *argv[]);


bool gQuietMode = false;
long gStrictLevel = 2;
bool gValidateInternals = true;
bool iniSettingsAreLegal(INIReader & iniReader, const mt::INIValueChecker & ivc, std::ostream &err);

int processContent(PublicNexusReader & nexusReader,
                   const char *,
                   std::ostream *os,
                   INIReader & currentAction);
MultiFormatReader * instantiateReader();

/* mtree INI checking */
bool iniSettingsAreLegal(INIReader & iniReader, const mt::INIValueChecker & ivc, std::ostream &err) {
    try{
        std::string value = iniReader.Get("action", "action", "LScore");
        ivc.parseActionAction(value); // triggers exception
        return true;
     } catch (mt::IllegalINIValueError & e) {
         err << "INI value error: " << e.what() << "\n";
         err.flush();
         return false;
     }
}
bool preDataINICheck(INIReader & iniReader, std::ostream & err) {
    mt::INIValueChecker ivc;
    return iniSettingsAreLegal(iniReader, ivc, err);
}
void INIReader::fill(mt::INIBasedSettings & ibs) {
    mt::INIValueChecker ivc;
    std::string value = this->Get("action", "action", "LScore");
    ibs.action = ivc.parseActionAction(value);
    value = this->Get("model", "ascertainment", "None");
    ibs.modelAsc = ivc.parseModelAscertainment(value);
}
/* end mtree INI checking */

namespace mt {




INIBasedSettings NCL2MT::configureBasedOnINI(MTInstance & , //mInstance,
                                               ::INIReader & iniReader,
                                               std::ostream & err) {
    INIBasedSettings ibs;
    iniReader.fill(ibs);
    err.flush();
    return ibs;
}

void NCL2MT::processTree(std::ostream *os,
            unsigned numTaxa,
            const NxsCharactersBlock * , //charsBlock,
            const NxsDiscreteDatatypeMapper * dataMapper,
            const NxsCDiscreteStateSet ** compressedMatrix,
            const vector<double> & patternWeights,
            const vector<int> &origToCompressed,
            const std::map<unsigned, std::set<unsigned> > & numStates2PatternIndexSet,
            const NxsSimpleTree & nxsTree,
            const ModelDescription & md,
            const INIBasedSettings & ibs) {
    assert(dataMapper != nullptr);
    vector<std::size_t> origToComp;
    for (auto otc : origToCompressed) {
        origToComp.push_back((std::size_t) otc);
    }
    std::map<unsigned, vector< vector<mt::char_state_t> > > rawPartMatrix;
    std::map<unsigned, mt::CharStateToPrimitiveInd > numStates2Cs2Pi;
    for (auto ns2pis : numStates2PatternIndexSet) {
        const unsigned numStates = ns2pis.first;
        const auto & patInds = ns2pis.second;
        vector< vector<mt::char_state_t> > & rawMatrix = rawPartMatrix[numStates];
        rawMatrix.resize(numTaxa);
        vector<mt::char_state_t> bogusChar;
        std::size_t currPartLen = patInds.size();
        if (md.GetAscBiasMode() == mt::ModelDescription::VAR_ONLY_NO_MISSING_ASC_BIAS) {
            currPartLen += numStates;
            for (auto i = 0U; i < numStates; ++i) {
                bogusChar.push_back(i);
            }
        }
        NxsCDiscreteStateSet maxStateCode = 0;
        for (auto i = 0U; i < numTaxa; ++i) {
            rawMatrix[i].reserve(currPartLen);
            for (auto j: patInds) {
                NxsCDiscreteStateSet r = compressedMatrix[i][j];
                //std::cerr << " compressedMatrix[" << i << "][" << j << "] << = " << (int)r << '\n';
                if (r < 0 || r > static_cast<NxsCDiscreteStateSet>(numStates)) {
                    r = static_cast<NxsCDiscreteStateSet>(numStates);
                }
                if (r > maxStateCode) {
                    maxStateCode = r;
                }
                rawMatrix[i].push_back((mt::char_state_t) r);
            }
            for (auto k : bogusChar) {
                rawMatrix[i].push_back(k);
            }
        }
        if (maxStateCode < (NxsCDiscreteStateSet) numStates) {
            maxStateCode = static_cast<NxsCDiscreteStateSet>(numStates);
        }
        unsigned numStateCodes = maxStateCode;
        auto & cs2pi = numStates2Cs2Pi[numStates];
        cs2pi.resize(numStateCodes + 1);
        std::set<mt::char_state_t> allFundStateCodes;
        for (auto i = 0U; i < numStateCodes; ++i) {
            vector<mt::char_state_t> v;
            for (auto xs : dataMapper->GetStateSetForCode(i)) {
                auto fsc = static_cast<mt::char_state_t>(xs);
                v.push_back(fsc);
                allFundStateCodes.insert(fsc);
            }
            cs2pi.SetStateCode(i, v);
        }
        vector<mt::char_state_t> allFundStateCodesVec(allFundStateCodes.begin(), allFundStateCodes.end());
        cs2pi.SetStateCode(numStateCodes, allFundStateCodesVec);
    }
    mt::MTInstance mtInstance(rawPartMatrix,
                              numStates2Cs2Pi,
                              patternWeights,
                              origToComp,
                              numStates2PatternIndexSet,
                              nxsTree,
                              md);
    doAnalysis(os, mtInstance, ibs);
}

} // namespace mt

int processContent(PublicNexusReader & nexusReader,
                   const char *gFilename,
                   std::ostream *os,
                   INIReader & iniReader) {
    BlockReaderList blocks = nexusReader.GetUsedBlocksInOrder();
    if (blocks.size() == 0) {
        cerr << "Error:\n No understandable content was found.\n";
        exit(1);
    }
    const unsigned numTaxaBlocks = nexusReader.GetNumTaxaBlocks();
    if (numTaxaBlocks != 1) {
        std::cerr << "Expecting a file with exactly 1 TAXA block, but found " << numTaxaBlocks << " in the file " << gFilename << ".\n";
        return 2;
    }
    NxsTaxaBlock * taxaBlock = nexusReader.GetTaxaBlock(0);
    const unsigned nCharBlocks = nexusReader.GetNumCharactersBlocks(taxaBlock);
    if (nCharBlocks != 1) {
        std::cerr << "Expecting a file with exactly 1 CHARACTERS/DATA block, but found " << nCharBlocks << " in the file " << gFilename << ".\n";
        return 3;
    }
    const unsigned nTreesBlocks = nexusReader.GetNumTreesBlocks(taxaBlock);
    if (nTreesBlocks != 1) {
        std::cerr << "Expecting a file with exactly 1 TREES block, but found " << nTreesBlocks << " in the file " << gFilename << ".\n";
        return 3;
    }
    const  NxsCharactersBlock * charBlock = nexusReader.GetCharactersBlock(taxaBlock, 0);
    assert(charBlock);
    std::vector<const NxsDiscreteDatatypeMapper *> mappers = charBlock->GetAllDatatypeMappers();
    if (mappers.size() != 1) {
        std::cerr << "Expecting an unmixed characters block, but found a matrix with datatype = mixed or a datatype with augmented symbols\n";
        return 4;
    }
    const NxsDiscreteDatatypeMapper * dm = mappers[0];
    ScopedTwoDMatrix<NxsCDiscreteStateSet> compressedMatrix;
    std::vector<unsigned> patternCounts;
    std::vector<double> patternWeights;
    bool hasWeights = true;
    std::vector<NxsCharacterPattern> compressedTransposedMatrix;
    std::vector<std::set<unsigned> > compressedIndexToOriginal;
    std::vector<int> originalIndexToCompressed;
    if (true) {
        if (true) {
            NxsCXXDiscreteMatrix cxxMat(*charBlock, false, 0L, false);
            hasWeights = cxxMat.hasWeights();
            NxsCompressDiscreteMatrix(cxxMat, compressedTransposedMatrix, &originalIndexToCompressed, &compressedIndexToOriginal);
        }
        std::vector<double> * wtsPtr = (hasWeights ? &patternWeights : 0L);
        NxsTransposeCompressedMatrix(compressedTransposedMatrix, compressedMatrix, &patternCounts, wtsPtr);
        patternWeights.clear();
        for (auto i : patternCounts) {
            patternWeights.push_back(double(i));
        }
    }
    const int MAX_NUM_STATES = 10; // TEMP. TODO, should be runtime var!
    std::map<unsigned, std::set<unsigned> > numStates2CharSet;
    for (int i = 2; i < MAX_NUM_STATES; ++i) {
        NxsUnsignedSet charIndexSet;
        auto x = NxsString::GetEscapedInt(i);
        x += " STATE CHARS";
        const auto nc = charBlock->GetIndexSet(x, &charIndexSet);
        if (nc > 0) {
            std::set<unsigned> patIndexSet;
            for (auto ci : charIndexSet) {
                auto pi = originalIndexToCompressed.at(ci);
                assert(pi >= 0);
                patIndexSet.insert(static_cast<unsigned>(pi));
            }
            std::cerr << nc << " chars and " << patIndexSet.size() << " patterns in \"" << x << "\"\n";
            numStates2CharSet[static_cast<unsigned>(i)] = patIndexSet;
        }
    }
    if (numStates2CharSet.empty()) {
        std::cerr << "Currently only works if the input file has charsets with names like 2_STATE_CHARS, 3_STATE_CHARS, etc.\n";
        assert(false);
        std::exit(1);
    }
    _DEBUG_VEC(patternWeights);
    NxsCDiscreteStateSet ** matrixAlias = compressedMatrix.GetAlias();
    const unsigned ntaxTotal =  charBlock->GetNTaxTotal();
    const  NxsTreesBlock * treesBlock = nexusReader.GetTreesBlock(taxaBlock, 0);
    mt::INIBasedSettings ibs;
    iniReader.fill(ibs);
    mt::ModelDescription md(ibs.modelAsc);
    mt::NCL2MT ncl2mt;
    for (unsigned nti = 0; nti < treesBlock->GetNumTrees(); ++nti) {
        const NxsSimpleTree nst(treesBlock->GetFullTreeDescription(nti), 1, 0.1, true);
        ncl2mt.processTree(os,
                           ntaxTotal,
                           charBlock,
                           dm,
                           (const NxsCDiscreteStateSet **) matrixAlias,
                           patternWeights,
                           originalIndexToCompressed,
                           numStates2CharSet,
                           nst,
                           md,
                           ibs);
    }
    return 0;
}


MultiFormatReader * instantiateReader() {
    MultiFormatReader * nexusReader = new MultiFormatReader(-1, NxsReader::WARNINGS_TO_STDERR);
    if (gQuietMode) {
        nexusReader->SetWarningOutputLevel(NxsReader::SKIPPING_CONTENT_WARNING);
    }
    if (gStrictLevel != 2) {
        nexusReader->SetWarningToErrorThreshold((int)NxsReader::FATAL_WARNING + 1 - (int) gStrictLevel);
    }
    NxsCharactersBlock * charsB = nexusReader->GetCharactersBlockTemplate();
    NxsDataBlock * dataB = nexusReader->GetDataBlockTemplate();
    charsB->SetAllowAugmentingOfSequenceSymbols(true);
    dataB->SetAllowAugmentingOfSequenceSymbols(true);
    NxsTreesBlock * treesB = nexusReader->GetTreesBlockTemplate();
    assert(treesB);
    if (gStrictLevel < 2) {
        treesB->SetAllowImplicitNames(true);
    }
    treesB->setValidateInternalNodeLabels(gValidateInternals);
    treesB->setAllowNumericInterpretationOfTaxLabels(true);
    if (gStrictLevel < 2) {
        NxsStoreTokensBlockReader *storerB =  nexusReader->GetUnknownBlockTemplate();
        assert(storerB);
        storerB->SetTolerateEOFInBlock(true);
    }
    nexusReader->conversionOutputRecord.addNumbersToDisambiguateNames = true;
    return nexusReader;
}

int processFilepath(
    const char * filename, // name of the file to be read
    std::ostream * os, // output stream to use (NULL for no output). Not that cerr is used to report errors.
    MultiFormatReader::DataFormatType fmt, // enum indicating the file format to expect.
    INIReader & iniReader) {
    assert(filename);
    int rc;
    try {
        MultiFormatReader * nexusReader;
        nexusReader = instantiateReader();
        if (!gQuietMode) {
            cerr << "Executing" << endl;
        }
        try {
            nexusReader->DemoteBlocks();
            nexusReader->ReadFilepath(filename, fmt);
            rc = processContent(*nexusReader, filename, os, iniReader);
        } catch(...) {
            nexusReader->DeleteBlocksFromFactories();
            delete nexusReader;
            throw;
        }
        nexusReader->DeleteBlocksFromFactories();
        delete nexusReader;
        return rc;
    } catch (const NxsException &x) {
        cerr << "Error:\n " << x.msg << endl;
        if (x.line > 0 || x.pos > 0) {
            cerr << "at line " << x.line << ", column (approximately) " << x.col << " (and file position "<< x.pos << ")" << endl;
        }
        return 2;
    }
}

int readFilepathAsNEXUS(const char *filename, MultiFormatReader::DataFormatType fmt, INIReader & iniReader) {
    if (!gQuietMode) {
        cerr << "[Reading " << filename << "     ]" << endl;
    }
    try {
        std::ostream * outStream = 0L;
        outStream = &cout;
        return processFilepath(filename, outStream, fmt, iniReader);
    } catch (...) {
        cerr << "Computing on " << filename << " failed (with an exception)" << endl;
        return 1;
    }
}

int readFilesListedIsFile(const char *masterFilepath, MultiFormatReader::DataFormatType fmt, INIReader & iniReader) {
    ifstream masterStream(masterFilepath, ios::binary);
    if (masterStream.bad()) {
        cerr << "Could not open " << masterFilepath << "." << endl;
        exit(3);
    }
    char filename[1024];
    while ((!masterStream.eof())  && masterStream.good()) {
        masterStream.getline(filename, 1024);
        if (strlen(filename) > 0 && filename[0] != '#') {
            int rc = readFilepathAsNEXUS(filename, fmt, iniReader);
            if (rc != 0) {
                return rc;
            }
        }
    }
    return 0;
}
extern const char * gExeName;
const char * gExeName = "mtree";
void printHelp(std::ostream & out) {
    out << "mtree takes reads a NEXUS file.\n";
    out << "\nThe most common usage is simply:\n    " << gExeName << " <path to NEXUS file>\n";
    out << "\nCommand-line flags:\n\n";
    out << "    -f<format> specifies the input file format expected:\n";
    out << "            -fnexus     NEXUS (this is also the default)\n";
    out << "            -faafasta   Amino acid data in fasta\n";
    out << "            -fdnafasta  DNA data in fasta\n";
    out << "            -frnafasta  RNA data in fasta\n";
    out << "        The complete list of format names that can follow the -f flag is:\n";
    std::vector<std::string> fmtNames =  MultiFormatReader::getFormatNames();
    for (std::vector<std::string>::const_iterator n = fmtNames.begin(); n != fmtNames.end(); ++n) {
        out << "            "<< *n << "\n";
    }
    out << "    -h help. on the command line shows this help message\n\n";
    out << "    -q quiet. suppress NCL status messages while reading files\n\n";
    out << "    -x do NOT validate internal labels in trees as taxa labels\n\n";
    out << "    -X do NOT treat numbers in trees as taxon numbers, treat them as arbitrary\n        labels (should not be used with NEXUS files).\n\n";
}

int do_main(int argc, char *argv[]) {
    NxsReader::setNCLCatchesSignals(true);
    MultiFormatReader::DataFormatType f(MultiFormatReader::NEXUS_FORMAT);
    std::string iniFilename;
    INIReader * iniReader = nullptr;
    for (int i = 1; i < argc; ++i) {
        const char * filepath = argv[i];
        const std::size_t slen = strlen(filepath);
        if (slen < 2 || filepath[0] != '-') {
            continue;
        }
        if (filepath[1] == 'h') {
            printHelp(cout);
            return 1;
        } else if (filepath[1] == 'q') {
            gQuietMode = true;
        } else if (filepath[1] == 'x') {
            gValidateInternals = false;
        } else if (filepath[1] == 'f') {
            f = MultiFormatReader::UNSUPPORTED_FORMAT;
            if (slen > 2) {
                std::string fmtName(filepath + 2, slen - 2);
                f =  MultiFormatReader::formatNameToCode(fmtName);
                if (f == MultiFormatReader::UNSUPPORTED_FORMAT) {
                    cerr << "Unknown format \"" << fmtName << "\" after -f\n" << endl;
                    printHelp(cerr);
                    return 3;
                }
            }
            if (f == MultiFormatReader::UNSUPPORTED_FORMAT) {
                cerr << "Expecting a format after -f\n" << endl;
                printHelp(cerr);
                return 2;
            }
        } else if (filepath[1] == 'm') {
            if (slen > 2) {
                if (!iniFilename.empty()) {
                    cerr << "Expecting one INI file with a -m flag.\n";
                    return 5;
                }
                iniFilename.assign(filepath + 2, slen - 2);
                iniReader = new INIReader(iniFilename.c_str());
                if (iniReader->ParseError() < 0) {
                    std::cerr << "Can't load \"" << iniFilename << "\"\n";
                    return 6;
                }
                if (!preDataINICheck(*iniReader, std::cerr)) {
                    std::cerr << "Exiting due to errors in \"" << iniFilename << "\"\n";
                    return 7;
                }
            } else {
                cerr << "Expecting an INI filepath after the -m flag\n";
                return 4;
            }
        }
    }
    if (iniReader == nullptr) {
        cerr << "Expecting an INI file specified with a -m flag\n";
        return 8;
    }
    bool readfile = false;
    for (int i = 1; i < argc; ++i) {
        const char * filepath = argv[i];
        const std::size_t slen = strlen(filepath);
        if (slen < 1)
            continue;
        if (strlen(filepath) > 2 && filepath[0] == '-' && filepath[1] == 'l') {
            readfile = true;
            int rc = readFilesListedIsFile(filepath+2, f, *iniReader);
            if (rc != 0) {
                return rc;
            }
        } else if (filepath[0] != '-') {
            readfile = true;
            int rc = readFilepathAsNEXUS(filepath, f, *iniReader);
            if (rc != 0) {
                return rc;
            }
        }
    }
    if (!readfile) {
        cerr << "Expecting the path to NEXUS file as the only command line argument!\n" << endl;
        printHelp(cerr);
        return 1;
    }
    return 0;
}

int main(int argc, char *argv[]) {
    int rc = do_main(argc, argv);
    return rc;
}

//
//  Based on normalizer.cpp from NCL (Nexus Class Library).
//
//  NCL is free software; you can redistribute it and/or modify
//  it under the terms of the GNU General Public License as published by
//  the Free Software Foundation; either version 2 of the License, or
//  (at your option) any later version.
//
//  NCL is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
//  GNU General Public License for more details.
//
//  You should have received a copy of the GNU General Public License
//  along with NCL; if not, write to the Free Software Foundation, Inc.,
//  59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
//
