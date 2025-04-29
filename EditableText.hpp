#pragma once

#include "Cairo.hpp"
#include <string>
#include <iostream>

START_NAMESPACE_DGL


class KeyboardEventHandler
{
public:
    class Callback
    {
    public:
        virtual ~Callback() {}
        virtual void keyPressed(SubWidget* widget, const Widget::KeyboardEvent& ev) = 0;
    };

    explicit KeyboardEventHandler(SubWidget* self);
    virtual ~KeyboardEventHandler();

    void setCallback(Callback* callback) noexcept;

    bool keyboardEvent(const Widget::KeyboardEvent& ev);

protected:
    virtual void stateChanged(const Widget::KeyboardEvent& ev);

private:
    struct PrivateData;
    PrivateData* const pData;

    DISTRHO_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR(KeyboardEventHandler)
};

struct KeyboardEventHandler::PrivateData {
    KeyboardEventHandler* const self;
    SubWidget* const widget;
    KeyboardEventHandler::Callback* userCallback;

    PrivateData(KeyboardEventHandler* const s, SubWidget* const w)
        : self(s),
          widget(w),
          userCallback(nullptr) {}

    bool handleKeyboardEvent(const Widget::KeyboardEvent& ev)
    {
        // Call the user's callback if set.
        if (userCallback != nullptr)
        {
            userCallback->keyPressed(widget, ev);
            return true;
        }

        // By default, let the self handle the event.
        self->stateChanged(ev);
        return false;
    }

    DISTRHO_DECLARE_NON_COPYABLE(PrivateData)
};

KeyboardEventHandler::KeyboardEventHandler(SubWidget* self)
    : pData(new PrivateData(this, self)) {}

KeyboardEventHandler::~KeyboardEventHandler()
{
    delete pData;
}

void KeyboardEventHandler::setCallback(Callback* const callback) noexcept
{
    pData->userCallback = callback;
}

bool KeyboardEventHandler::keyboardEvent(const Widget::KeyboardEvent& ev)
{
    return pData->handleKeyboardEvent(ev);
}

void KeyboardEventHandler::stateChanged(const Widget::KeyboardEvent&)
{
    // Default implementation does nothing.
    // Subclasses can override this to handle state changes.
}


class EditableText : public CairoSubWidget,
                            public ButtonEventHandler,
                            public ButtonEventHandler::Callback,
                            public KeyboardEventHandler
{
public:
    class Callback
    {
    public:
        virtual ~Callback() {}
        virtual void setTF(std::string newText) = 0;
        virtual void focus(EditableText* current) = 0;
        virtual void setGain(std::string gainStr) = 0;
        virtual void generateIfTracking(void) = 0;
    };

    explicit EditableText(SubWidget* const parent, std::string type)
        : CairoSubWidget(parent),
          ButtonEventHandler(this),
          KeyboardEventHandler(this),
          fCallback(nullptr),
          fText(""),
          fIsFocused(false),
          fCoeffCount(1),
          fType(type)
    {
        for (int i = 0; i < 100; i++) {
            fPoint[i] = false;
            fMinus[i] = false;
        }
        ButtonEventHandler::setCallback(this);
    }

    explicit EditableText(TopLevelWidget* const parent, std::string type)
        : CairoSubWidget(parent),
          ButtonEventHandler(this),
          KeyboardEventHandler(this),
          fCallback(nullptr),
          fText(""),
          fIsFocused(false),
          fCoeffCount(1),
          fType(type)
    {
        for (int i = 0; i < 100; i++) {
            fPoint[i] = false;
            fMinus[i] = false;
        }
        ButtonEventHandler::setCallback(this);
    }

    void setText(const std::string& text)
    {
        fText = text;
        repaint();
    }

    const std::string& getText() const
    {
        return fText;
    }

    void setFocused(bool newState) {
        fIsFocused = newState;
    }

    bool isFocused() {
        return fIsFocused;
    }

    void setCallback(Callback* const callback) noexcept
    {
        fCallback = callback;
    }

protected:
    void onCairoDisplay(const CairoGraphicsContext& context) override
    {
        cairo_t* cr = context.handle;

        // Draw background
        cairo_set_source_rgb(cr, 1.0, 1.0, 1.0);
        cairo_rectangle(cr, 0, 0, getWidth(), getHeight());
        cairo_fill(cr);

        // Draw border (highlight if focused)
        if (this->isFocused())
            cairo_set_source_rgb(cr, 0.0, 0.0, 1.0); // Blue border for focus
        else
            cairo_set_source_rgb(cr, 0.0, 0.0, 0.0); // Black border otherwise

        cairo_rectangle(cr, 0, 0, getWidth(), getHeight());
        cairo_stroke(cr);

        // Draw text
        cairo_set_source_rgb(cr, 0.0, 0.0, 0.0);
        cairo_select_font_face(cr, "Sans", CAIRO_FONT_SLANT_NORMAL, CAIRO_FONT_WEIGHT_NORMAL);
        cairo_set_font_size(cr, 14);

        cairo_move_to(cr, 5, getHeight() - 10); // Adjust text position
        cairo_show_text(cr, fText.c_str());
    }

    bool onMouse(const MouseEvent& ev) override
    {
        std::cout << "onMouse: " << this << "\n" << std::flush;
        if (fCallback)
            fCallback->focus(nullptr);
        return ButtonEventHandler::mouseEvent(ev);
    }

    void buttonClicked(SubWidget* widget, int button) override
    {
        std::cout << "buttonClicked: " << this << "\n" << std::flush;
        fIsFocused = true;
        repaint();
        if (fCallback)
            fCallback->focus(this);
    }

    bool onKeyboard(const KeyboardEvent& ev) override
    {
        if (this->isFocused() && ev.press)
        {
            if (fType == "G") {
                if ((ev.key >= '0' && ev.key <= '9') || (ev.key == '.' && !fPoint[fCoeffCount - 1])) {
                    if (ev.key == '.') fPoint[fCoeffCount - 1] = true;
                    fText += ev.key;
                    repaint();
                }
                else if (ev.key == 8 && !fText.empty()) // Handle backspace
                {
                    if (fText[fText.length() - 1] == '.') fPoint[fCoeffCount - 1] = false;
                    if (fText[fText.length() - 1] == ',') fCoeffCount--;
                    if (fText[fText.length() - 1] == '-') fMinus[fCoeffCount - 1] = false;
                    fText.pop_back();
                    repaint();
                }
                else if (ev.key == 13) // Use enter key to apply gain
                {
                    this->setFocused(false);
                    fCallback->setGain(fText);
                    fText = "";
                    fPoint[fCoeffCount - 1] = false;
                    repaint();
                }
            } else {
                if ((ev.key >= '0' && ev.key <= '9') || (ev.key == '.' && !fPoint[fCoeffCount - 1]) || (ev.key == '-' && !fMinus[fCoeffCount - 1] && (!(std::strcmp(fText.c_str(), "")) || fText[fText.length() - 1] == ',')) || (ev.key == ',')) // Allow digits and decimal point
                {
                    if (ev.key == '.') fPoint[fCoeffCount - 1] = true;
                    if (ev.key == ',') {
                        fPoint[++fCoeffCount - 1] = false;
                        fMinus[++fCoeffCount - 1] = false;
                    }
                    if (ev.key == '-') fMinus[fCoeffCount - 1] = true;
                    fText += ev.key;
                    repaint();
                }
                else if (ev.key == 8 && !fText.empty()) // Handle backspace
                {
                    if (fText[fText.length() - 1] == '.') fPoint[fCoeffCount - 1] = false;
                    if (fText[fText.length() - 1] == ',') fCoeffCount--;
                    if (fText[fText.length() - 1] == '-') fMinus[fCoeffCount - 1] = false;
                    fText.pop_back();
                    repaint();
                }
                else if (ev.key == 13) // Enter key
                {
                    this->setFocused(false); // Lose focus on Enter
                    repaint();
                }
                if (fCallback) {
                    if (fType != " ")
                        fCallback->setTF(fType + fText);
                    else
                        fCallback->generateIfTracking();
                }
            }
        }
        return KeyboardEventHandler::keyboardEvent(ev);
    }

private:
    Callback* fCallback;
    std::string fText;
    bool fIsFocused;
    bool fPoint[100];
    bool fMinus[100];
    int fCoeffCount;
    std::string fType;
};

END_NAMESPACE_DGL
