/*
 * DISTRHO Plugin Framework (DPF)
 * Copyright (C) 2012-2021 Filipe Coelho <falktx@falktx.com>
 * Copyright (C) 2019-2021 Jean Pierre Cimalando <jp-dev@inbox.ru>
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

#include "DistrhoUI.hpp"

#include "Artwork.hpp"
#include "DemoWidgetBanner.hpp"
#include "DemoWidgetClickable.hpp"
#include  "EditableText.hpp"
#include <cstdio>
#include <iostream>

START_NAMESPACE_DISTRHO

// We need a few classes from DGL.
using DGL_NAMESPACE::CairoGraphicsContext;
using DGL_NAMESPACE::CairoImage;
using DGL_NAMESPACE::CairoImageKnob;
using DGL_NAMESPACE::CairoImageSwitch;

// And from ourselves
using DGL_NAMESPACE::DemoWidgetBanner;
using DGL_NAMESPACE::EditableText;

class FildesUI : public UI,
                       public CairoImageKnob::Callback,
                       public CairoImageSwitch::Callback,
                       public DemoWidgetClickable::Callback,
                       public EditableText::Callback
{
    ScopedPointer<CairoImageKnob> fKnob;
    ScopedPointer<CairoImageSwitch> fButton;
    ScopedPointer<DemoWidgetClickable> fWidgetClickable;
    ScopedPointer<EditableText> fCoeffT, fCoeffB;
    ScopedPointer<DGL::SubWidget> fContainer;

public:
    FildesUI()
    {
        freopen("output.txt","w",stdout);
        std::cout << "start\n" << std::flush;
        CairoImage knobSkin;
        knobSkin.loadFromPNG(Artwork::knobData, Artwork::knobDataSize);

        fWidgetClickable = new DemoWidgetClickable(this);
        fWidgetClickable->setAbsolutePos(100, 100);
        fWidgetClickable->setSize(50, 50);
        fWidgetClickable->setCallback(this);
        fWidgetClickable->setId(kParameterTriState);
        fWidgetClickable->setVisible(false);

        fKnob = new CairoImageKnob(this, knobSkin);
        fKnob->setAbsolutePos(10, 100);
        fKnob->setSize(80, 80);
        fKnob->setCallback(this);
        fKnob->setId(kParameterKnob);
        fKnob->setVisible(false);

        CairoImage buttonOn, buttonOff;
        buttonOn.loadFromPNG(Artwork::buttonOnData, Artwork::buttonOnDataSize);
        buttonOff.loadFromPNG(Artwork::buttonOffData, Artwork::buttonOffDataSize);

        fButton = new CairoImageSwitch(this, buttonOff, buttonOn);
        fButton->setAbsolutePos(100, 160);
        fButton->setSize(60, 35);
        fButton->setCallback(this);
        fButton->setId(kParameterButton);
        fButton->setVisible(false);

        // Initialize in constructor
        fCoeffT = new EditableText(this, true);
        fCoeffT->setAbsolutePos(50, 50);
        fCoeffT->setSize(500, 30);
        fCoeffT->setCallback(this);

        fCoeffB = new EditableText(this, false);
        fCoeffB->setAbsolutePos(50, 100);
        fCoeffB->setSize(500, 30);
        fCoeffB->setCallback(this);

        // we can use this if/when our resources are scalable, for now they are PNGs
        const double scaleFactor = getScaleFactor();
        setSize(800,400);
        if (scaleFactor != 1.0)
            setSize(400 * scaleFactor, 800 * scaleFactor);
    }

protected:
    void onCairoDisplay(const CairoGraphicsContext& context)
    {
        cairo_t* const cr = context.handle;
        cairo_set_source_rgb(cr, 0.541, 0.576, 1);
        cairo_paint(cr);

        cairo_set_source_rgb(cr, 0.0, 0.0, 0.0); // Set text color
        cairo_select_font_face(cr, "Sans", CAIRO_FONT_SLANT_NORMAL, CAIRO_FONT_WEIGHT_NORMAL);
        cairo_set_font_size(cr, 20);

        cairo_move_to(cr, 10, 30);
        cairo_show_text(cr, "Transfer Function:"); // Static text

    }

    // we can use this if/when our resources are scalable, for now they are PNGs
    void onResize(const ResizeEvent& ev) override
    {
        UI::onResize(ev);

        const double scaleFactor = getScaleFactor();

        fWidgetClickable->setSize(50*scaleFactor, 50*scaleFactor);
        fWidgetClickable->setAbsolutePos(100*scaleFactor, 100*scaleFactor);

        fKnob->setSize(80*scaleFactor, 80*scaleFactor);
        fKnob->setAbsolutePos(10*scaleFactor, 100*scaleFactor);

        fButton->setSize(60*scaleFactor, 35*scaleFactor);
        fButton->setAbsolutePos(100*scaleFactor, 160*scaleFactor);
    }

    void parameterChanged(const uint32_t index, const float value) override
    {
        switch (index)
        {
        case kParameterKnob:
            fKnob->setValue(value);
            break;
        case kParameterTriState:
            fWidgetClickable->setColorId(static_cast<int>(value + 0.5f));
            break;
        case kParameterButton:
            fButton->setDown(value > 0.5f);
            break;
        }
    }
    
    void focus(EditableText* current) override {
        if (current == fCoeffT) {
            fCoeffB->setFocused(false);
        } else if (current == fCoeffB) {
            fCoeffT->setFocused(false);
        } else {
            fCoeffB->setFocused(false);
            fCoeffT->setFocused(false);
        }
        repaint();
    }

    void setTF(std::string newText) override {
        std::cout << "setTF\n" << std::flush;
        if (newText[0] == 'T') setState("topTF", &newText[1]);
        else if (newText[0] == 'B') setState("bottomTF", &newText[1]);
    }

    void stateChanged(const char*, const char*) {
        return;
    }

    void demoWidgetClicked(DemoWidgetClickable*, const uint8_t colorId) override
    {
        setParameterValue(kParameterTriState, colorId);
    }

    void imageKnobDragStarted(CairoImageKnob*) override
    {
        editParameter(kParameterKnob, true);
    }

    void imageKnobDragFinished(CairoImageKnob*) override
    {
        editParameter(kParameterKnob, false);
    }

    void imageKnobValueChanged(CairoImageKnob*, const float value) override
    {
        setParameterValue(kParameterKnob, value);
    }

    void imageSwitchClicked(CairoImageSwitch*, bool down) override
    {
        setParameterValue(kParameterButton, down ? 1.f : 0.f);
    }
};

UI* createUI()
{
    return new FildesUI;
}

END_NAMESPACE_DISTRHO
