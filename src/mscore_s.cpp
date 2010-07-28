/*
 * Copyright (c) 2003-2006 Fred Hutchinson Cancer Research Center
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

#include "stdafx.h"
#include "msequence.h"
#include "mspectrum.h"
#include "xmlparameter.h"
#include "mscore_s.h"

// Factory instance, registers itself with the mscoremanager.
static mscorefactory_s factory;
    
mscorefactory_s::mscorefactory_s()
{
    mscoremanager::register_factory("s-score", this);
}

mplugin* mscorefactory_s::create_plugin()
{
    return new mscore_s();
}

mscore_s::mscore_s(void)
{
    m_dScale = 0.05;    
}

mscore_s::~mscore_s(void)
{
}

bool mscore_s::clear()
{
    m_vmiType.clear();
    return true;
}

/*
 * allows score object to issue warnings, or set variable based on xml.
 */
bool mscore_s::load_param(XmlParameter &_x)
{
    if (!mscore::load_param(_x))
        return false;

    string strValue;
    string strKey = "s-score, histogram scale";
    if(_x.get(strKey,strValue))    {
        m_dScale = atof(strValue.c_str());
    }

    return true;
}

void mscore_s::prescore(const size_t _i)
{
    mscore::prescore(_i);

    m_used.clear();
}
/*
 * add_mi does the work necessary to set up an mspectrum object for modeling. 
 *   - an entry in the m_State object is made for the parent ion M+H
 * once an mspectrum has been added, the original mspectrum is no longer
 * needed for modeling, as all of the work associated with a spectrum
 * is only done once, prior to modeling sequences.
 */
bool mscore_s::add_mi(mspectrum &_s)
{
    if (!mscore::add_mi(_s))
        return false;

    vmiType vType;
    MIType uType;
    double fTotal = 0.0;

    vector<mi>::iterator itMI = _s.m_vMI.begin();
    vector<mi>::iterator itEnd = _s.m_vMI.end();
    while (itMI != itEnd) {
        uType.m_lM = (long) (itMI->m_fM + 0.5);
        uType.m_fI = log(itMI->m_fI);
        vType.push_back(uType);

        itMI++;
    }

    m_vmiType.push_back(vType);
    return true;
}
/*
 * mconvert converts from mass and charge to integer ion value
 * for mi vector.
 */
unsigned long mscore_s::mconvert(double _m, const long _c)
{
    const double fZ = (double)_c;

    double dMass = (fZ*m_pSeqUtilFrag->m_dProton + _m)/fZ;
    return (long)(dMass + 0.5);
}
/*
 * sfactor returns a factor applied to the final convolution score.
 */
double mscore_s::sfactor()
{
    return 10.0 / sqrt(strlen(m_pSeq));
}
/*
 * report_score formats a hyper score for output.
 */
void mscore_s::report_score(char* buffer, float hyper)
{
    sprintf(buffer, "%d",(int) (hyper + 0.5));
}

/*
 * dot is the fundamental logic for scoring a peptide with a mass spectrum.
 * the mass spectrum is determined by the value of m_lId, which is its index
 * number in the m_vsmapMI vector. the sequence is represented by the values
 * that are currently held in m_plSeq (integer masses).
 */
double mscore_s::dot(unsigned long *_v)
{
    float fScore = 0.0;
    unsigned long lCount = 0;

    vector<MIType>::iterator itType = m_vmiType[m_lId].begin();
    vector<MIType>::const_iterator itEnd = m_vmiType[m_lId].end();
    for (unsigned long a = 0; m_plSeq[a] != 0; a++) {

        while(itType != itEnd && itType->m_lM < m_plSeq[a])
            itType++;

        if (itType != itEnd && itType->m_lM == m_plSeq[a]) {
            if (itType->m_fI > 0 && m_used.find(itType->m_lM) == m_used.end())    {
                lCount++;
                m_used.insert(itType->m_lM);
                fScore += itType->m_fI;
            }
        }
    }

    *_v = lCount;
    return (fScore);
}
