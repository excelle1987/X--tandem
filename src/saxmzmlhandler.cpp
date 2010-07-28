/************************************************************
 *              SAXMzmlHandler.cpp
 * Adapted from SAXMzdataHandler.cpp
 * August 2008
 * Ronald Beavis
 *
 * Premiere version janvier 2005
 * Patrick Lacasse
 * placasse@mat.ulaval.ca
 *
 * 3/11/2005 (Brendan MacLean): Use eXpat SAX parser, and create SAXSpectraHandler
 *
 * November 2005
 * Fredrik Levander 
 * Fredrik.Levander@elmat.lth.se 
 * A few changes to handle MzData 1.05.
 *
 * Updated to handle version 1.04 and 1.05. (Rob Craig)
 *
 *
 * See http://psidev.sourceforge.net/ms/#mzdata for
 * mzData schema information.
 *
 * Inspired by DtaSAX2Handler.cpp
 * copyright            : (C) 2002 by Pedrioli Patrick, ISB, Proteomics
 * email                : ppatrick@systemsbiology.org
 * Artistic License granted 3/11/2005
 *******************************************************/

#include "stdafx.h"
#include "saxmzmlhandler.h"

SAXMzmlHandler::SAXMzmlHandler( vector<mspectrum>& _vS, mspectrumcondition& _sC, mscore& _m)
	: SAXSpectraHandler(_vS, _sC, _m)
{
	m_bInmzArrayBinary = false;
	m_bInintenArrayBinary = false;
	m_bInData = false;
	m_bInMsLevel2 = false;
	m_bNetworkData = false;
	m_bLowPrecision = false;
	m_bGaml = false;
}

SAXMzmlHandler::~SAXMzmlHandler()
{
}

void SAXMzmlHandler::startElement(const XML_Char *el, const XML_Char **attr)
{
	// Let's see what	element	just started
	// We	use:
	//   spectrum	for	the	id number
	//   acqInstrument for the msLevel
	//   cvParam for the polarity	and	the	mz of the parent ion
	//   data	for	the	values of the spectrum

	if (isElement("spectrum", el))
	{
		reset();	// Clean up	for	the	next scan
		m_scanNum =	atoi(getAttrValue("index", attr));
		m_tId = m_scanNum;
		while(m_sId.find(m_tId) != m_sId.end())	{
			m_tId++;
		}
		m_sId.insert(m_tId);
		int iLevel = atoi(getAttrValue("msLevel", attr));
		if(iLevel == 2)	{
			m_bInMsLevel2 = true;
		}
		m_peaksCount = atoi(getAttrValue("defaultArrayLength", attr));
	}
	else if (isElement("softwareParam", el))	{
		const char* name = getAttrValue("name", attr);
		const char* accession = getAttrValue("accession", attr);
		const char* version = getAttrValue("version", attr);
	}
	else if (isElement("cvParam", el))
	{
		const char* name = getAttrValue("name", attr);
		const char* accession = getAttrValue("accession", attr);
		const char* value = getAttrValue("value", attr);
		if ((!strcmp(name, "ms level") || !strcmp(accession,"MS:1000511")) && !strcmp(value,"2")){
			m_bInMsLevel2 = true;
		}
		else if(!strcmp(name, "charge state") || !strcmp(accession,"MS:1000041"))	{
			m_precursorCharge = atoi(value);
		}
		else if(!strcmp(name, "m/z") || !strcmp(accession,"MS:1000040"))	{
			m_precursorMz =	atof(value);
		}
		else if(!strcmp(name, "64-bit float") || !strcmp(accession,"MS:1000523"))	{
			m_bLowPrecision = false;
		}
		else if(!strcmp(name, "32-bit float") || !strcmp(accession,"MS:1000521"))	{
			m_bLowPrecision = true;
		}
		else if(!strcmp(name, "m/z array") || !strcmp(accession,"MS:1000514"))	{
			m_bInData = true;
			m_bInmzArrayBinary = true;
			m_bInintenArrayBinary = false;
		}
		else if(!strcmp(name, "intensity array") || !strcmp(accession,"MS:1000515"))	{
			m_bInData = true;
			m_bInintenArrayBinary = true;
			m_bInmzArrayBinary = false;
		}
	}
	if(isElement("binary", el))	{
		m_strData.clear();
		m_bInData = true;
	}
	if(isElement("binaryDataArray", el))	{
		m_strData.clear();
		if(atoi(getAttrValue("arrayLength", attr)) > 0)	{
			m_peaksCount = atoi(getAttrValue("arrayLength", attr));
		}
		m_bInData = true;
	}
}


void SAXMzmlHandler::endElement(const XML_Char *el)
{
	if(isElement("binary", el))	{
		processData();
		m_bInintenArrayBinary = false;
		m_bInmzArrayBinary = false;
		m_bInData = false;
	}
	else if(isElement("spectrum", el) && m_bInMsLevel2)
	{
		pushSpectrum();
		m_bInMsLevel2 = false;
	}
}

void SAXMzmlHandler::characters(const XML_Char *s, int len)
{
	if(m_bInmzArrayBinary && m_bInMsLevel2 && m_bInData){
		m_strData.append(s, len);
	}
	if(m_bInintenArrayBinary && m_bInMsLevel2 && m_bInData){
		m_strData.append(s, len);
	}
}

void SAXMzmlHandler::processData()
{
	if(m_bInmzArrayBinary && m_bInMsLevel2 && m_bInData)
	{
		pushPeaks(m_bInmzArrayBinary, m_bInintenArrayBinary);
	}
	else if(m_bInintenArrayBinary && m_bInMsLevel2 && m_bInData)
	{
		pushPeaks(m_bInmzArrayBinary, m_bInintenArrayBinary);
	}

	m_strData.clear();
}
