#include <cassert>
#include <map>
#include "ncl/ncl.h"
#include "ncl/nxsblock.h"
#include "ncl/nxspublicblocks.h"
#include "ncl/nxscxxdiscretematrix.h"
#include "ncl/nxsmultiformat.h"
#include "mt_tree.h"
using namespace std;

void ncl2mt(unsigned numTaxa,
              const NxsDiscreteDatatypeMapper * dataMapper,
             const NxsCDiscreteStateSet ** compressedMatrix,
             const vector<double> & patternWeights,
             const vector<int> &origToCompressed,
             const NxsSimpleTree & tree,
             const mt::ModelDescription & md);


bool gQuietMode = false;
long gStrictLevel = 2;
bool gValidateInternals = true;
enum ProcessActionsEnum {
    SCORE_ACTION=0
};


int processContent(PublicNexusReader & nexusReader, const char *, std::ostream *os, ProcessActionsEnum currentAction);
MultiFormatReader * instantiateReader();

void ncl2mt(unsigned numTaxa,
            const NxsDiscreteDatatypeMapper * dataMapper,
             const NxsCDiscreteStateSet ** compressedMatrix,
             const vector<double> & patternWeights,
             const vector<int> &origToCompressed,
             const NxsSimpleTree & nxsTree,
             const mt::ModelDescription & md) {
    assert(dataMapper != nullptr);
    const unsigned numRealChars = patternWeights.size();
    unsigned firstPartLength = patternWeights.size();
    const unsigned numStates = dataMapper->GetNumStates();
    vector<unsigned> origToComp;
    for (auto otc : origToCompressed) {
        origToComp.push_back((unsigned) otc);
    }
    vector<mt::char_state_t> bogusChar;
    if (md.GetAscBiasMode() == mt::ModelDescription::VAR_ONLY_NO_MISSING_ASC_BIAS) {
        firstPartLength += numStates;
        for (auto i = 0U; i < numStates; ++i) {
            bogusChar.push_back(i);
        }
    }
    vector<vector<mt::char_state_t> > rawMatrix(numTaxa);
    NxsCDiscreteStateSet maxStateCode = 0;
    for (auto i = 0U; i < numTaxa; ++i) {
        rawMatrix[i].reserve(firstPartLength);
        for (auto j = 0U; j < numRealChars; ++j) {
            NxsCDiscreteStateSet r = compressedMatrix[i][j];
            if (r < 0) {
                r = numStates;
            } else if (r > maxStateCode) {
                maxStateCode = r;
            }
            rawMatrix[i].push_back((mt::char_state_t) compressedMatrix[i][j]);
        }
        for (auto k : bogusChar) {
            rawMatrix[i].push_back(k);
        }
    }
    vector<mt::char_state_t *> rowPtrs(numTaxa);
    for (auto i = 0U; i < numTaxa; ++i) {
        rowPtrs[i] = &(rawMatrix[i][0]);
    }
    if (maxStateCode < (NxsCDiscreteStateSet) numStates) {
        maxStateCode = numStates;
    }
    unsigned nsc = maxStateCode;
    mt::CharStateToPrimitiveInd cs2pi(nsc);
    for (auto i = 0U; i < nsc; ++i) {
        vector<mt::char_state_t> v;
        for (auto xs : dataMapper->GetStateSetForCode(i)) {
            v.push_back(static_cast<mt::char_state_t>(xs));
        }
        cs2pi.SetStateCode(i, v);
    }

    std::vector<unsigned> partLengths(1, firstPartLength);
    mt::PartitionedMatrix partMat(numTaxa, partLengths, origToComp);
    partMat.fillPartition(0, const_cast<const mt::char_state_t**>(&(rowPtrs[0])), &cs2pi);
    unsigned numNodes = 2 * numTaxa - 1;
    mt::Tree tree(numNodes, numTaxa);
    std::map<const NxsSimpleNode *, unsigned> ncl2nodeNumber;
    std::vector<const NxsSimpleNode *> pre = nxsTree.GetPreorderTraversal();
    unsigned internalIndex = numTaxa;
    for (std::vector<const NxsSimpleNode *>::iterator ndIt = pre.begin(); ndIt != pre.end(); ++ndIt) {
        const NxsSimpleNode *nd = *ndIt;
        std::cout << "address = " << (long) nd << " is leaf = " << nd->IsTip() << " index = " << nd->GetTaxonIndex() << " parent address = " << (long) nd->GetEdgeToParent().GetParent() << std::endl;
        unsigned num;
        if (nd->IsTip()) {
            num = nd->GetTaxonIndex();
        } else {
            num = internalIndex++;
        }
        mt::Node * treeNode = tree.GetNode(num);
        const NxsSimpleNode * par = nd->GetEdgeToParent().GetParent();
        if (par == nullptr) {
            tree.SetRoot(treeNode);
        } else {
            assert (ncl2nodeNumber.find(par) ! = ncl2nodeNumber.end());
            unsigned parNodeNumber = ncl2nodeNumber[par];
            mt::Node * parNode = tree.GetNode(parNodeNumber);
            parNode->AddChild(treeNode, nd->GetEdgeToParent().GetDblEdgeLen());
        }
        ncl2nodeNumber[nd] = num;
    }
    unsigned numRateCats = 1;
    for (auto li = 0U; li < numTaxa; ++li) {
        mt::Node * leaf = tree.GetLeaf(li);
        assert(leaf);
        for (auto j = 0U; j < partMat.GetNumPartitions(); ++j) {
            leaf->SetData(j, (void *) partMat.GetLeafCharacters(j, li));
            leaf->SetWork(j, (void *) new mt::LeafWork(maxStateCode + 1, numStates, numRateCats));
        }
    }
    for (auto li = numTaxa; li < numNodes; ++ li) {
        mt::Node * nd = tree.GetNode(li);
        for (auto j = 0U; j < partMat.GetNumPartitions(); ++j) {
            nd->SetWork(j, (void *) new mt::InternalNodeWork(firstPartLength, numStates, numRateCats));
        }
    }
    mt::CharModel * cm;
    if (md.GetAscBiasMode() == mt::ModelDescription::VAR_ONLY_NO_MISSING_ASC_BIAS) {
        cm = new mt::MkVarNoMissingAscCharModel(numStates);
    } else {
        cm = new mt::MkCharModel(numStates);
    }
    try {
        doAnalysis(tree, *cm);
    } catch (...) {
        delete cm;
        throw;
    }
    delete cm;
}

int processContent(PublicNexusReader & nexusReader,
                   const char *gFilename,
                   std::ostream *os,
                   ProcessActionsEnum ) {
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
    bool hasIntWeights = true;
    std::vector<NxsCharacterPattern> compressedTransposedMatrix;
    std::vector<std::set<unsigned> > compressedIndexToOriginal;
    std::vector<int> originalIndexToCompressed;
    if (true) {
        if (true) {
            NxsCXXDiscreteMatrix cxxMat(*charBlock, false, 0L, false);
            hasWeights = cxxMat.hasWeights();
            hasIntWeights = cxxMat.hasIntWeights();
            NxsCompressDiscreteMatrix(cxxMat, compressedTransposedMatrix, &originalIndexToCompressed, &compressedIndexToOriginal);
        }
       std::vector<double> * wtsPtr = (hasWeights ? &patternWeights : 0L);
       NxsTransposeCompressedMatrix(compressedTransposedMatrix, compressedMatrix, &patternCounts, wtsPtr);
    }
    NxsCDiscreteStateSet ** matrixAlias = compressedMatrix.GetAlias();
    const unsigned ntaxTotal =  charBlock->GetNTaxTotal();
    const unsigned numPatterns = patternCounts.size();
    *os << "#NEXUS\nBEGIN DATA;\n\tDimensions ntax = " << ntaxTotal << " nchar = " << numPatterns << ";\n\t";
    charBlock->WriteFormatCommand(*os);
    *os << "Matrix\n";
    const unsigned width = taxaBlock->GetMaxTaxonLabelLength();
    for (unsigned i = 0; i < ntaxTotal; i++) {
        const std::string currTaxonLabel = NxsString::GetEscaped(taxaBlock->GetTaxonLabel(i));
        *os << currTaxonLabel;
        unsigned currTaxonLabelLen = (unsigned)currTaxonLabel.size();
        unsigned diff = width - currTaxonLabelLen;
        for (unsigned k = 0; k < diff+5; k++) {
            *os << ' ';
        }
        NxsCDiscreteStateSet * matrixRow = matrixAlias[i];
        for (unsigned j = 0; j < numPatterns; ++j) {
            *os << (int) matrixRow[j] << ' ';
        }
        *os << '\n';
    }
    const char * sp = (hasWeights ? " " : " * ");
    unsigned ind = 0;
    *os << ";\nEND;\nBEGIN ASSUMPTIONS;\n\tWTSET" << sp << " counts ( vector ) =";
    for (std::vector<unsigned>::const_iterator cIt = patternCounts.begin(); cIt != patternCounts.end(); ++cIt, ++ind) {
        *os << ' ' << *cIt;
    }
    *os << ";\n";
    if (hasWeights) {
        *os << "\tWTSET * sum_of_weights ( vector ) =";
        if (hasIntWeights) {
            for (std::vector<double>::const_iterator cIt = patternWeights.begin(); cIt != patternWeights.end(); ++cIt, ++ind) {
                int w = int(0.01 + *cIt);
                *os << ' ' << w;
            }
        } else {
            for (std::vector<double>::const_iterator cIt = patternWeights.begin(); cIt != patternWeights.end(); ++cIt, ++ind) {
                *os << ' ' << *cIt;
            }
        }
        *os << ";\n";
    }
    *os << "END;\n";
    const  NxsTreesBlock * treesBlock = nexusReader.GetTreesBlock(taxaBlock, 0);
    mt::ModelDescription md(mt::ModelDescription::VAR_ONLY_NO_MISSING_ASC_BIAS); //@TODO should be run-time setting
    //mt::ModelDescription md(mt::ModelDescription::NO_ASC_BIAS); //@TODO should be run-time setting
    for (unsigned nti = 0; nti < treesBlock->GetNumTrees(); ++nti) {
        const NxsSimpleTree nst(treesBlock->GetFullTreeDescription(nti), 1, 0.1, true);
        ncl2mt(ntaxTotal, dm, (const NxsCDiscreteStateSet **) matrixAlias, patternWeights, originalIndexToCompressed, nst, md);
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
    ProcessActionsEnum currentAction) {
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
            rc = processContent(*nexusReader, filename, os, currentAction);
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

int readFilepathAsNEXUS(const char *filename, MultiFormatReader::DataFormatType fmt, ProcessActionsEnum currentAction) {
    if (!gQuietMode) {
        cerr << "[Reading " << filename << "     ]" << endl;
    }
    try {
        std::ostream * outStream = 0L;
        outStream = &cout;
        return processFilepath(filename, outStream, fmt, currentAction);
    } catch (...) {
        cerr << "Normalizing of " << filename << " failed (with an exception)" << endl;
        return 1;
    }
}

int readFilesListedIsFile(const char *masterFilepath, MultiFormatReader::DataFormatType fmt, ProcessActionsEnum currentAction) {
    ifstream masterStream(masterFilepath, ios::binary);
    if (masterStream.bad()) {
        cerr << "Could not open " << masterFilepath << "." << endl;
        exit(3);
    }
    char filename[1024];
    while ((!masterStream.eof())  && masterStream.good()) {
        masterStream.getline(filename, 1024);
        if (strlen(filename) > 0 && filename[0] != '#') {
            int rc = readFilepathAsNEXUS(filename, fmt, currentAction);
            if (rc != 0) {
                return rc;
            }
        }
    }
    return 0;
}

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
    ProcessActionsEnum currentAction= SCORE_ACTION;
    NxsReader::setNCLCatchesSignals(true);
    MultiFormatReader::DataFormatType f(MultiFormatReader::NEXUS_FORMAT);
    for (int i = 1; i < argc; ++i) {
        const char * filepath = argv[i];
        const unsigned slen = strlen(filepath);
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
        }
    }
    bool readfile = false;
    for (int i = 1; i < argc; ++i) {
        const char * filepath = argv[i];
        const unsigned slen = strlen(filepath);
        if (slen < 1)
            continue;
        if (strlen(filepath) > 2 && filepath[0] == '-' && filepath[1] == 'l') {
            readfile = true;
            int rc = readFilesListedIsFile(filepath+2, f, currentAction);
            if (rc != 0) {
                return rc;
            }
        } else if (filepath[0] != '-') {
            readfile = true;
            int rc = readFilepathAsNEXUS(filepath, f, currentAction);
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
