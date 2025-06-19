#include "DistrhoPlugin.hpp"

#include <string.h>
#include <iostream>
#include <cstdio>
#include <bitset>
#include <atomic>
#include <mutex>

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

    }

    const char* getLabel() const override
    {
        return "fildes";
    }

    const char* getDescription() const override
    {
        return "An educational filter designer";
    }

    const char* getMaker() const override
    {
        return "Declan McMahon";
    }

    const char* getLicense() const override
    {
        return "MIT";
    }

    uint32_t getVersion() const override
    {
        return d_version(1, 0, 0);
    }

    void initAudioPort(const bool input, const uint32_t index, AudioPort& port) override
    {
        // treat meter audio ports as stereo
        port.groupId = kPortGroupMono;

        // everything else is as default
        Plugin::initAudioPort(input, index, port);
    }

    void initParameter(const uint32_t index, Parameter& parameter) override
    {
       
    }

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

        char inputCopy[1024];
        strncpy(inputCopy, str, sizeof(inputCopy) - 1);
        inputCopy[sizeof(inputCopy) - 1] = '\0';

        char* token = strtok(inputCopy, ",");
        while (token != NULL && index < 100) {
            double value = (token[0] == '\0') ? 0.0f : strtof(token, NULL);
            outputArray[index] = value;

            if (value != 0.0f) {
                *numCoeffs = (int)index;
            }

            token = strtok(NULL, ",");
            ++index;
        }

        for (; index < 100; ++index) {
            outputArray[index] = 0.0f;
        }

        fPrevChangei = fRuni;
        fTFChanged = true;
    }

    void setState(const char* key, const char* value) override {
        if (!std::strcmp(key, "topTF")) {
            processTFString(value, true);
        } else if (!std::strcmp(key, "bottomTF")) {
            processTFString(value, false);
        }
    }

    void setParameterValue(const uint32_t index, const float value) override
    {
        
    }

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

                uint32_t bits;
                std::memcpy(&bits, &outputs[n][i], sizeof(bits));

                if (bits == 0xFFC00000) {
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
