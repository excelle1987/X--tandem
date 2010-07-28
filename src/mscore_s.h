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

#include "mscore.h"
#include <set>

/*
 * mscore_s implements a very simple scoring algorithm.
 */
class mscore_s : public mscore
{
protected:
    friend class mscorefactory_s;

    mscore_s(void);    // Should only be created through mscorefactory_s

public:
    virtual ~mscore_s(void);

public:
    virtual bool load_param(XmlParameter &_x); // allows score object to issue warnings,
                                               // or set variables based on xml 
    virtual void prescore(const size_t _i); // called before scoring

    virtual bool add_mi(mspectrum &_s);

    virtual double sfactor(); // factor applied to final convolution score
    virtual unsigned long mconvert(double _m, const long _c); // convert mass to integer ion m/z for mi vector
    virtual void report_score(char* _buff, float _h); // format hyper score for output

    virtual bool clear();
    /* Thu, Apr 22, 2010 at 9:26 AM on spctools-discuss, Brendan M said:
    I took a quick look at the code, and this function appears to be
    related to new handling for phosphorylation.  It looked to me like you
    could just implement the function in k-score, and return zero without
    any ill effects, unless you are trying to use the new phosphorylation
    functionality.  Then all bets are off, without someone really
    understanding the meaning of this function and how the multipliers
    applied to the score when ion_check is non-zero apply to k-score, and
    then testing a real implementation.
    But, I think you can compile and use it with just a zero returning
    function.
    */
    float ion_check(long unsigned int foo, size_t bar) {
        return 0.0f;  // cross your fingers
    }


protected:
    virtual double dot(unsigned long *_v); // this is where the real scoring happens

protected:
    vector<vmiType> m_vmiType;

    set<unsigned long> m_used;
};

/*
 * mscorefactory_s implements a factory for creating mscore_s instances.
 */
class mscorefactory_s : public mpluginfactory
{
public:
    mscorefactory_s();

    virtual mplugin* create_plugin();
};
