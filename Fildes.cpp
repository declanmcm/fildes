/*
 * DISTRHO Plugin Framework (DPF)
 * Copyright (C) 2019-2021 Jean Pierre Cimalando <jp-dev@inbox.ru>
 * Copyright (C) 2012-2024 Filipe Coelho <falktx@falktx.com>
 *
 * Permission to use, copy, modify, and/or distribute this software for any purpose with
 * or without fee is hereby granted, provided that the above copyright notice and this
 * permission notice appear in all copies.
 *
 * THE SOFTWARE IS PROVIDED "AS IS" AND THE AUTHOR DISCLAIMS ALL WARRANTIES WITH REGARD
 * TO THIS SOFTWARE INCLUDING ALL IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS. IN
 * NO EVENT SHALL THE AUTHOR BE LIABLE FOR ANY SPECIAL, DIRECT, INDIRECT, OR CONSEQUENTIAL
 * DAMAGES OR ANY DAMAGES WHATSOEVER RESULTING FROM LOSS OF USE, DATA OR PROFITS, WHETHER
 * IN AN ACTION OF CONTRACT, NEGLIGENCE OR OTHER TORTIOUS ACTION, ARISING OUT OF OR IN
 * CONNECTION WITH THE USE OR PERFORMANCE OF THIS SOFTWARE.
 */

#include "DistrhoPlugin.hpp"

#include <string.h>
#include <iostream>
#include <cstdio>

START_NAMESPACE_DISTRHO

// HP IIR
//double acoeff[]={-0.000012029930817286773,0.00028060999630655074,-0.0030974027000735135,0.02152174301835309,-0.1050673162817373,0.3832808442810133,-1.0747462874237625,2.373876325566058,-4.122317208384148,5.69304181789048,-6.026599020577933,4.928492002867855,-2.683565830256412,1};
//double bcoeff[]={1,13,78,286,715,1287,1716,1716,1287,715,286,78,13,1};

// LP IIR
//double acoeff[]={0.04946174058699353,0.12611476775484237,0.4025931735415687,0.2006010749469854,1.1660419715456853,-0.16666928763817113,1};
//double bcoeff[]={1,-6,15,-20,15,-6,1};

// LP FIR
/*double acoeff[] = {
    0.000000000000000000,
    -0.000005622290253718,
    -0.000015952613718074,
    -0.000019164857121407,
    0.000000000000000000,
    0.000061151997198953,
    0.000189005632131738,
    0.000413777450584962,
    0.000771052143228499,
    0.001300980553345611,
    0.002046732813480528,
    0.003052205201402669,
    0.004359063037605834,
    0.006003285535379580,
    0.008011453006166186,
    0.010397073738067305,
    0.013157280105781339,
    0.016270226109782456,
    0.019693489417184757,
    0.023363721007878532,
    0.027197698740929570,
    0.031094834376003570,
    0.034941065856108781,
    0.038613948413921587,
    0.041988650142707258,
    0.044944470273701805,
    0.047371440083511518,
    0.049176543178285452,
    0.050289106802944697,
    0.050664968287519985,
    0.050289106802944704,
    0.049176543178285452,
    0.047371440083511518,
    0.044944470273701798,
    0.041988650142707251,
    0.038613948413921600,
    0.034941065856108795,
    0.031094834376003567,
    0.027197698740929577,
    0.023363721007878549,
    0.019693489417184764,
    0.016270226109782463,
    0.013157280105781339,
    0.010397073738067312,
    0.008011453006166198,
    0.006003285535379582,
    0.004359063037605840,
    0.003052205201402673,
    0.002046732813480528,
    0.001300980553345611,
    0.000771052143228498,
    0.000413777450584962,
    0.000189005632131737,
    0.000061151997198953,
    0.000000000000000000,
    -0.000019164857121407,
    -0.000015952613718074,
    -0.000005622290253718,
    0.000000000000000000
};*/

// HP FIR
double acoeff[] = {
    0.000000000000000000,
    -0.000025695329810175,
    -0.000109001586731170,
    -0.000261383834469497,
    -0.000497284187006939,
    -0.000834034055775751,
    -0.001291444409523049,
    -0.001891070645768172,
    -0.002655173095064218,
    -0.003605418278370960,
    -0.004761388110941483,
    -0.006138982684280024,
    -0.007748815659053334,
    -0.009594708572398906,
    -0.011672390798343837,
    -0.013968505240990876,
    -0.016460006299700354,
    -0.019114016903342568,
    -0.021888186575993940,
    -0.024731564041054279,
    -0.027585967536853457,
    -0.030387805705078649,
    -0.033070273559075526,
    -0.035565823484471772,
    -0.037808792095928072,
    -0.039738051371860414,
    -0.041299547699896848,
    -0.042448595699213505,
    -0.043151804843812071,
    0.956611464609619988,
    -0.043151804843812071,
    -0.042448595699213505,
    -0.041299547699896848,
    -0.039738051371860407,
    -0.037808792095928072,
    -0.035565823484471785,
    -0.033070273559075540,
    -0.030387805705078645,
    -0.027585967536853460,
    -0.024731564041054296,
    -0.021888186575993951,
    -0.019114016903342579,
    -0.016460006299700354,
    -0.013968505240990888,
    -0.011672390798343854,
    -0.009594708572398910,
    -0.007748815659053343,
    -0.006138982684280031,
    -0.004761388110941481,
    -0.003605418278370961,
    -0.002655173095064216,
    -0.001891070645768171,
    -0.001291444409523047,
    -0.000834034055775751,
    -0.000497284187006940,
    -0.000261383834469497,
    -0.000109001586731170,
    -0.000025695329810175,
    0.000000000000000000
};

class Fildes : public Plugin
{
    double topTF[100];
    double bottomTF[100];
    double outMem[100];
	double inMem[100];
	double gain;
    int numCoeffsTop, numCoeffsBottom, outi, ini, xz1, xz2, yz1, yz2;

public:
    Fildes()
        : Plugin(kParameterCount, 0, 0)
    {
        std::memset(topTF, 0, sizeof(topTF));
        std::memset(bottomTF, 0, sizeof(bottomTF));
        bottomTF[0] = 1;
        std::memset(outMem, 0, sizeof(outMem));
        std::memset(inMem, 0, sizeof(inMem));
        numCoeffsTop = numCoeffsBottom = ini = outi = xz1 = xz2 = yz1 = yz2 = 0;
		std::cout << "Fildes\n" << std::flush;

		//std::memcpy(topTF, acoeff, sizeof(acoeff));
		gain = 1;
		//std::memcpy(bottomTF, bcoeff, sizeof(bcoeff));

    }

   /**
      Get the plugin label.@n
      This label is a short restricted name consisting of only _, a-z, A-Z and 0-9 characters.
    */
    const char* getLabel() const override
    {
        return "cairo_ui";
    }

   /**
      Get an extensive comment/description about the plugin.
    */
    const char* getDescription() const override
    {
        return "Cairo DPF Example";
    }

   /**
      Get the plugin author/maker.
    */
    const char* getMaker() const override
    {
        return "DISTRHO";
    }

   /**
      Get the plugin license (a single line of text or a URL).@n
      For commercial plugins this should return some short copyright information.
    */
    const char* getLicense() const override
    {
        return "ISC";
    }

   /**
      Get the plugin version, in hexadecimal.
    */
    uint32_t getVersion() const override
    {
        return d_version(1, 0, 0);
    }

   /**
      Initialize the audio port @a index.@n
      This function will be called once, shortly after the plugin is created.
    */
    void initAudioPort(const bool input, const uint32_t index, AudioPort& port) override
    {
        // treat meter audio ports as stereo
        port.groupId = kPortGroupMono;

        // everything else is as default
        Plugin::initAudioPort(input, index, port);
    }

   /**
      Initialize the parameter @a index.@n
      This function will be called once, shortly after the plugin is created.
    */
    void initParameter(const uint32_t index, Parameter& parameter) override
    {
       
    }

   /**
      Get the current value of a parameter.@n
      The host may call this function from any context, including realtime processing.
    */
    float getParameterValue(const uint32_t index) const override
    {
        return 0;
    }

    void processTFString(const char* str, bool type) {

        size_t index = 0;
        
        int *numCoeffs;
        double *outputArray;
        if (type) {
            numCoeffs = &numCoeffsTop;
            outputArray = topTF;
        }
        else {
            numCoeffs = &numCoeffsBottom;
            outputArray = &bottomTF[1];
        }
        *numCoeffs = -1;

        // Create a copy of the input string to tokenize
        char inputCopy[1024]; // Ensure this is large enough for your input strings
        strncpy(inputCopy, str, sizeof(inputCopy) - 1);
        inputCopy[sizeof(inputCopy) - 1] = '\0'; // Ensure null-termination

        // Tokenize the string using commas as the delimiter
        char* token = strtok(inputCopy, ",");
        while (token != NULL && index < 100) {
            // Convert the token to a float, treating empty tokens as 0
            double value = (token[0] == '\0') ? 0.0f : strtof(token, NULL);
            outputArray[index] = value;

            // Update the highest non-zero index
            if (value != 0.0f) {
                *numCoeffs = (int)index;
            }

            // Move to the next token
            token = strtok(NULL, ",");
            ++index;
        }

        // Fill remaining elements with 0 if the array is larger than the input
        for (; index < 100; ++index) {
            outputArray[index] = 0.0f;
        }

        for (int i = 0; i < 4; i++) {
            std::cout << topTF[i] << "\n" << std::flush;
            std::cout << bottomTF[i] << "\n" << std::flush;
        }
    }

    void setState(const char* key, const char* value) override {
        std::cout << "setState\n" << std::flush;
        if (!std::strcmp(key, "topTF")) {
            processTFString(value, true);
        } else if (!std::strcmp(key, "bottomTF")) {
            processTFString(value, false);
        }
    }

   /**
      Change a parameter value.@n
      The host may call this function from any context, including realtime processing.@n
      When a parameter is marked as automatable, you must ensure no non-realtime operations are performed.
      @note This function will only be called for parameter inputs.
    */
    void setParameterValue(const uint32_t index, const float value) override
    {
        
    }

   /**
      Run/process function for plugins without MIDI input.
      @note Some parameters might be null if there are no audio inputs or outputs.
    */
    void run(const float** const inputs, float** const outputs, const uint32_t frames) override
    {

        if (inputs == nullptr || outputs == nullptr || inputs[0] == nullptr || outputs[0] == nullptr) {
            return;
        }

        int found = 0;
        double outputCopy[frames];

        for (int i = 0; i < frames; i++) {
            double res = 0;
            
            if (i > 98) {
                for (int j = 0; j < 100; j++) {
                    res += (double)inputs[0][i - j] * topTF[j];
                    if (j != 0) res -= (double)outputCopy[i - j] * bottomTF[j];
                }
                if (frames - i <= 100) {
                    inMem[(++ini) % 100] = (double)inputs[0][i];
                    outMem[(++outi) % 100] = res;
                }
            } else {
                for (int j = 0; j < i + 1; j++) {
                    res += (double)inputs[0][i - j] * topTF[j];
                    if (j != 0) res -= (double)outputCopy[i - j] * bottomTF[j];
                }
                for (int j = 0; j < 99 - i; j++) {
                    int inIndex = (ini - j) % 100;
                    int outIndex = (outi - j) % 100;
                    res += inMem[inIndex] * topTF[i + 1 + j];
                    if (i + 1 + j != 0) res -= outMem[outIndex] * bottomTF[i + 1 + j];
                }
            }

            outputs[0][i] = res * gain;
            outputCopy[i] = res * gain;
        }
    }

};

Plugin* createPlugin()
{
    return new Fildes();
}

END_NAMESPACE_DISTRHO
