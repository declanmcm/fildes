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
#include <bitset>

START_NAMESPACE_DISTRHO

class Fildes : public Plugin
{
    double topTF[100];
    double bottomTF[100];
    double fNewTop[100];
    double fNewBottom[100];
    double fOldTop[100];
    double fOldBottom[100];
    double outMem[2][100];
	double inMem[2][100];
	double gain;
    double alpha;
    int numCoeffsTop, numCoeffsBottom, outi, ini, xz1, xz2, yz1, yz2, fRuni, fPrevChangei, fInterpolIter, fMaxValCount;
    bool fTFChanged;

public:
    Fildes()
        : Plugin(kParameterCount, 0, 0)
    {
        std::memset(topTF, 0, sizeof(topTF));
        std::memset(bottomTF, 0, sizeof(bottomTF));
        std::memset(fNewTop, 0, sizeof(fNewTop));
        std::memset(fNewBottom, 0, sizeof(fNewBottom));
        std::memset(fOldTop, 0, sizeof(fOldTop));
        std::memset(fOldBottom, 0, sizeof(fOldBottom));
        bottomTF[0] = 1;
        topTF[0] = 1;
        fNewBottom[0] = 1;
        std::memset(outMem, 0, sizeof(outMem));
        std::memset(inMem, 0, sizeof(inMem));
        numCoeffsTop = numCoeffsBottom = ini = outi = fMaxValCount = xz1 = xz2 = yz1 = yz2 = 0;
		std::cout << "Fildes\n" << std::flush;
        fTFChanged = false;

        fRuni = 0;
        fPrevChangei = 0;
        alpha = 0.01;

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
        return "fildes";
    }

   /**
      Get an extensive comment/description about the plugin.
    */
    const char* getDescription() const override
    {
        return "An educational filter designer";
    }

   /**
      Get the plugin author/maker.
    */
    const char* getMaker() const override
    {
        return "Declan McMahon";
    }

   /**
      Get the plugin license (a single line of text or a URL).@n
      For commercial plugins this should return some short copyright information.
    */
    const char* getLicense() const override
    {
        return "MIT";
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
            outputArray = fNewTop;
        }
        else {
            numCoeffs = &numCoeffsBottom;
            outputArray = &fNewBottom[1];
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

        // for (int i = 0; i < 100; i++) {
        //     fOldTop[i] = topTF[i];
        //     fOldBottom[i] = bottomTF[i];
        // }

        for (int i = 0; i < 4; i++) {
            std::cout << topTF[i] << "\n" << std::flush;
            std::cout << bottomTF[i] << "\n" << std::flush;
        }

        fPrevChangei = fRuni;
        fTFChanged = true;
    }

    void setState(const char* key, const char* value) override {
        //std::cout << "setState\n" << std::flush;
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

        for (int n = 0; n < 2; n++) {
            for (int i = 0; i < frames; i++) {
                double res = 1e-20;

                fRuni++;
                if(fRuni > 99999) {
                    fPrevChangei -= fRuni;
                    fRuni = 0;
                }

                if (fTFChanged && (fRuni - fPrevChangei > 100)) {
                    fPrevChangei = fRuni;

                    std::memcpy(topTF, fNewTop, sizeof(fNewTop));
                    std::memcpy(bottomTF, fNewBottom, sizeof(fNewBottom));
                    std::memset(inMem, 0, sizeof(inMem));
                    std::memset(outMem, 0, sizeof(outMem));
                    fTFChanged = false;
                }
                
                if (i > 98) {
                    for (int j = 0; j < 100; j++) {
                        res += (double)inputs[n][i - j] * topTF[j];
                        if (j != 0) res -= (double)outputCopy[i - j] * bottomTF[j];
                    }
                    if (frames - i <= 100) {
                        inMem[n][(++ini) % 100] = (double)inputs[n][i];
                        outMem[n][(++outi) % 100] = res;
                    }
                } else {
                    for (int j = 0; j < i + 1; j++) {
                        res += (double)inputs[n][i - j] * topTF[j];
                        if (j != 0) res -= (double)outputCopy[i - j] * bottomTF[j];
                    }
                    for (int j = 0; j < 99 - i; j++) {
                        int inIndex = (ini - j) % 100;
                        int outIndex = (outi - j) % 100;
                        res += inMem[n][inIndex] * topTF[i + 1 + j];
                        if (i + 1 + j != 0) res -= outMem[n][outIndex] * bottomTF[i + 1 + j];
                    }
                }

                outputs[n][i] = static_cast<float>(res);
                outputCopy[i] = res;

                //std::cerr << outputs[n][i] << std::endl;

                uint32_t bits;
                std::memcpy(&bits, &outputs[n][i], sizeof(bits));

                if (bits == 0xFFC00000) {
                    //std::cout << "TOO LOUD" << std::endl;
                    fMaxValCount++;
                    if (fMaxValCount > 100) {
                        fMaxValCount = 0;
                        std::memset(inMem, 0, sizeof(inMem));
                        std::memset(outMem, 0, sizeof(outMem));
                    }
                }
            }
        }
    }

};

Plugin* createPlugin()
{
    return new Fildes();
}

END_NAMESPACE_DISTRHO
