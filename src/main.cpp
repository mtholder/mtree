#include "ncl/ncl.h"
#include "ncl/nxsblock.h"
#include "ncl/nxspublicblocks.h"
#include "ncl/nxsmultiformat.h"
#include <cassert>
using namespace std;
//#include "ncl/nxscxxdiscretematrix.h"

bool gQuietMode = false;
std::ofstream gCommonFileStream;
std::ostream * gCommonOstream = 0L;
long gStrictLevel = 2;
bool gUnderscoresToSpaces = false;
bool gValidateInternals = true;
bool gTreesViaInMemoryStruct = true;
long gInterleaveLen = -1;
bool blocksReadInValidation = false;
bool gSuppressingNameTranslationFile = false;
bool gAllowNumericInterpretationOfTaxLabels = true;
TranslatingConventions gTranslatingConventions;

void processContent(PublicNexusReader & nexusReader, ostream *os, ProcessActionsEnum currentAction);
MultiFormatReader * instantiateReader();
MultiFormatReader * gNexusReader = NULL;


////////////////////////////////////////////////////////////////////////////////
// Takes NxsReader that has successfully read a file, and processes the
//	information stored in the reader.
//
// The caller is responsibel for calling DeleteBlocksFromFactories() to clean
//	up (if the reader uses the factory API).
////////////////////////////////////////////////////////////////////////////////
void processContent(PublicNexusReader & nexusReader, ostream *os, ProcessActionsEnum currentAction) {
	BlockReaderList blocks = nexusReader.GetUsedBlocksInOrder();
	if (blocks.size() == 0) {
		cerr << "Error:\n No understandable content was found.\n";
		exit(1);
	}
}


MultiFormatReader * instantiateReader()
{
	MultiFormatReader * nexusReader = new MultiFormatReader(-1, NxsReader::WARNINGS_TO_STDERR);
	if (gQuietMode)
		nexusReader->SetWarningOutputLevel(NxsReader::SKIPPING_CONTENT_WARNING);
	if (gStrictLevel != 2)
		nexusReader->SetWarningToErrorThreshold((int)NxsReader::FATAL_WARNING + 1 - (int) gStrictLevel);
	if (gUnderscoresToSpaces) 
		nexusReader->SetCoerceUnderscoresToSpaces(true);
	NxsCharactersBlock * charsB = nexusReader->GetCharactersBlockTemplate();
	NxsDataBlock * dataB = nexusReader->GetDataBlockTemplate();
	charsB->SetAllowAugmentingOfSequenceSymbols(true);
	dataB->SetAllowAugmentingOfSequenceSymbols(true);
	if (gInterleaveLen > 0)
		{
		assert(charsB);
		charsB->SetWriteInterleaveLen(gInterleaveLen);
		dataB->SetWriteInterleaveLen(gInterleaveLen);
		}

	NxsTreesBlock * treesB = nexusReader->GetTreesBlockTemplate();
	assert(treesB);
	if (gStrictLevel < 2)
		treesB->SetAllowImplicitNames(true);
	treesB->SetWriteFromNodeEdgeDataStructure(gTreesViaInMemoryStruct);
	treesB->setValidateInternalNodeLabels(gValidateInternals);
	if (!gValidateInternals) {
		gTranslatingConventions.treatNodeLabelsAsStrings = true;
	}
	treesB->setAllowNumericInterpretationOfTaxLabels(gAllowNumericInterpretationOfTaxLabels);
	if (gAltNexus)
		treesB->setWriteTranslateTable(false);
	if (gStrictLevel < 2)
		{
		NxsStoreTokensBlockReader *storerB =  nexusReader->GetUnknownBlockTemplate();
		assert(storerB);
		storerB->SetTolerateEOFInBlock(true);
		}
	nexusReader->conversionOutputRecord.addNumbersToDisambiguateNames = true;
	
	if (gSuppressingNameTranslationFile)
		nexusReader->conversionOutputRecord.writeNameTranslationFile = false;
	return nexusReader;
}


int processFilepath(
	const char * filename, // name of the file to be read
	ostream * , // output stream to use (NULL for no output). Not that cerr is used to report errors.
	MultiFormatReader::DataFormatType fmt, // enum indicating the file format to expect.
	ProcessActionsEnum ) // enum that is passed on to processContent to indicate what should be done with the content of the file.
	{
	assert(filename);
	try
		{
		MultiFormatReader * nexusReader;
		nexusReader = instantiateReader();

		if (!gQuietMode)
			cerr << "Executing" << endl;
		try {
			nexusReader->DemoteBlocks();
			nexusReader->ReadFilepath(filename, fmt);
			processContent(*nexusReader, os, currentAction);
			}
		catch(...)
			{
			nexusReader->DeleteBlocksFromFactories();
			delete nexusReader;
			throw;
			}
		nexusReader->DeleteBlocksFromFactories();
		delete nexusReader;
		return 0;
		}
	catch (const NxsException &x)
		{
		cerr << "Error:\n " << x.msg << endl;
		if (x.line > 0 || x.pos > 0)
			cerr << "at line " << x.line << ", column (approximately) " << x.col << " (and file position "<< x.pos << ")" << endl;
		return 2;
		}
	}

/*! \returns 0 on success*/
int readFilepathAsNEXUS(const char *filename, MultiFormatReader::DataFormatType fmt, ProcessActionsEnum currentAction)
	{
	if (!gQuietMode)
		cerr << "[Reading " << filename << "	 ]" << endl;
	try {
		ostream * outStream = 0L;
		if (currentAction != VALIDATE_ONLY)
			outStream = &cout;
		return processFilepath(filename, outStream, fmt, currentAction);

		}
	catch (...)
		{
		cerr << "Normalizing of " << filename << " failed (with an exception)" << endl;
		return 1;
		}
	}

/*! \returns 0 on success*/
int readFilesListedIsFile(const char *masterFilepath, MultiFormatReader::DataFormatType fmt, ProcessActionsEnum currentAction)
	{
	ifstream masterStream(masterFilepath, ios::binary);
	if (masterStream.bad())
		{
		cerr << "Could not open " << masterFilepath << "." << endl;
		exit(3);
		}
	char filename[1024];
	while ((!masterStream.eof())  && masterStream.good())
		{
		masterStream.getline(filename, 1024);
		if (strlen(filename) > 0 && filename[0] != '#')
			{
			int rc = readFilepathAsNEXUS(filename, fmt, currentAction);
			if (rc != 0)
				return rc;
			}
		}
	return 0;
	}

const char * gExeName = "mtree";

void printHelp(ostream & out) {
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
	for (std::vector<std::string>::const_iterator n = fmtNames.begin(); n != fmtNames.end(); ++n)
		{
		out << "            "<< *n << "\n";
		}
	out << "    -h help. on the command line shows this help message\n\n";
	out << "    -q quiet. suppress NCL status messages while reading files\n\n";
	out << "    -x do NOT validate internal labels in trees as taxa labels\n\n";
	out << "    -X do NOT treat numbers in trees as taxon numbers, treat them as arbitrary\n        labels (should not be used with NEXUS files).\n\n";
	}

int do_main(int argc, char *argv[])
	{
	NxsReader::setNCLCatchesSignals(true);
	MultiFormatReader::DataFormatType f(MultiFormatReader::NEXUS_FORMAT);
	for (int i = 1; i < argc; ++i)
		{
		const char * filepath = argv[i];
		const unsigned slen = strlen(filepath);
		if (slen < 2 || filepath[0] != '-')
			continue;
		if (filepath[1] == 'h')
			{
			printHelp(cout);
			return 1;
			}
		else if (filepath[1] == 'q')
			gQuietMode = true;
		else if (filepath[1] == 'x')
			gValidateInternals = false;
		else if (filepath[1] == 'X')
			gAllowNumericInterpretationOfTaxLabels = false;
		else if (filepath[1] == 'u')
			gUnderscoresToSpaces = true;
		else if (filepath[1] == 'j')
			gSuppressingNameTranslationFile = true;
		else if (filepath[1] == 'f')
			{
			f = MultiFormatReader::UNSUPPORTED_FORMAT;
			if (slen > 2)
				{
				std::string fmtName(filepath + 2, slen - 2);
				f =  MultiFormatReader::formatNameToCode(fmtName);
				if (f == MultiFormatReader::UNSUPPORTED_FORMAT)
					{
					cerr << "Unknown format \"" << fmtName << "\" after -f\n" << endl;
					printHelp(cerr);
					return 3;
					}
				}
			if (f == MultiFormatReader::UNSUPPORTED_FORMAT)
				{
				cerr << "Expecting a format after -f\n" << endl;
				printHelp(cerr);
				return 2;
				}
			}
		}
	bool readfile = false;
	for (int i = 1; i < argc; ++i)
		{
		const char * filepath = argv[i];
		const unsigned slen = strlen(filepath);
		if (slen < 1)
			continue;
		if (strlen(filepath) > 2 && filepath[0] == '-' && filepath[1] == 'l')
			{
			readfile = true;
			int rc = readFilesListedIsFile(filepath+2, f, currentAction);
			if (rc != 0)
				return rc;
			}
		else if (filepath[0] != '-')
			{
			readfile = true;
			int rc = readFilepathAsNEXUS(filepath, f, currentAction);
			if (rc != 0)
				return rc;
			}
		}

	if (!readfile)
		{
		cerr << "Expecting the path to NEXUS file as the only command line argument!\n" << endl;
		printHelp(cerr);
		return 1;
		}
	return 0;
}

int main(int argc, char *argv[]) {
	int rc = do_main(argc, argv);
	if (gCommonOstream != 0L && gCommonOstream == &gCommonFileStream)
		gCommonFileStream.close();
	return rc;
}

//
//	Based on normalizer.cpp from NCL (Nexus Class Library).
//
//	NCL is free software; you can redistribute it and/or modify
//	it under the terms of the GNU General Public License as published by
//	the Free Software Foundation; either version 2 of the License, or
//	(at your option) any later version.
//
//	NCL is distributed in the hope that it will be useful,
//	but WITHOUT ANY WARRANTY; without even the implied warranty of
//	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
//	GNU General Public License for more details.
//
//	You should have received a copy of the GNU General Public License
//	along with NCL; if not, write to the Free Software Foundation, Inc.,
//	59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
//
