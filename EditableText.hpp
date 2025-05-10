#pragma once

#include "Cairo.hpp"
#include <string>
#include <iostream>
#include <chrono>
#include <string>
#include <thread>
 #include <atomic>
#include <mutex>

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

    int cursorPosition;
    std::vector<double> charXPositions;
    std::string fText, fPreviousText;
    int fCursorX, fCursorY;
    bool fIsLayoutValid;
    cairo_t* fCr;
    // Debouncing for text changes
    std::atomic<bool> fTextChanged;
    std::chrono::time_point<std::chrono::steady_clock> fLastTextChangeTime;
    const int fTextChangeThrottleMs = 50; // Minimum ms between text change callbacks

public:
    class Callback
    {
    public:
        virtual ~Callback() {}
        virtual void setTF(std::string newText) = 0;
        virtual void focus(EditableText* current) = 0;
        virtual void setGain(std::string gainStr) = 0;
        virtual void setMaxA(std::string maxAStr) = 0;
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
        cursorPosition = 0;
        for (int i = 0; i < 100; i++) {
            fPoint[i] = false;
            fMinus[i] = false;
        }
        fIsLayoutValid = false;
        fLastTextChangeTime = std::chrono::steady_clock::now();
        ButtonEventHandler::setCallback(this);
        
        startTextChangeThread();
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
        fLastTextChangeTime = std::chrono::steady_clock::now();
        ButtonEventHandler::setCallback(this);
        
        // Start the text change handler thread
        startTextChangeThread();
    }

    ~EditableText()
    {
        // Signal the thread to stop
        fTextChangeThreadRunning = false;
        if (fTextChangeThread.joinable()) {
            fTextChangeThread.join();
        }
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
        fCr = cr;

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

        drawText(cr, 5, getHeight() - 10);
    }

    void drawText(cairo_t* cr, double x, double y) {
        // Store the beginning position
        charXPositions.clear();
        charXPositions.push_back(x);

        fPreviousText = fText;
        
        // Select font before measuring
        cairo_select_font_face(cr, "monospace", CAIRO_FONT_SLANT_NORMAL, CAIRO_FONT_WEIGHT_NORMAL);
        cairo_set_font_size(cr, 14);
        
        // Calculate position for each character
        double currentX = x;
        for (size_t i = 0; i < fText.length(); i++) {
            cairo_text_extents_t extents;
            // Measure single character
            std::string charStr = fText.substr(i, 1);
            cairo_text_extents(cr, charStr.c_str(), &extents);
            
            // Draw character
            cairo_move_to(cr, currentX, y);
            cairo_show_text(cr, charStr.c_str());
            
            // Move to next position and store it
            currentX += extents.x_advance;
            charXPositions.push_back(currentX);
        }
        
        // Draw cursor if needed
        if (cursorPosition >= 0 && cursorPosition <= fText.length() && this->isFocused()) {
            cairo_set_source_rgb(cr, 0, 0, 0);
            cairo_move_to(cr, charXPositions[cursorPosition], y - 14); // Adjust height based on font
            cairo_line_to(cr, charXPositions[cursorPosition], y + 2);
            cairo_stroke(cr);
        }

        fIsLayoutValid = true;
    }

    int getCharIndexAtPosition(double x, double y) {
        // Find which character boundary the click is closest to
        for (size_t i = 0; i < charXPositions.size() - 1; i++) {
            if (x >= charXPositions[i] && x < charXPositions[i+1]) {
                // Check if click is closer to the left or right character boundary
                double leftDist = x - charXPositions[i];
                double rightDist = charXPositions[i+1] - x;
                return (leftDist < rightDist) ? i : i + 1;
            }
        }
        
        // If beyond the end of text, return the last position
        if (x >= charXPositions.back()) {
            return fText.length();
        }
        
        return -1;
    }

    bool onMouse(const MouseEvent& ev) override
    {
        if (fCallback)
            fCallback->focus(nullptr);

        if (ev.button == 1) {
            if (ev.press) {
                int newCursorPos = getCharIndexAtPosition(ev.pos.getX(), ev.pos.getY());
                if (newCursorPos >= 0)
                    cursorPosition = newCursorPos;
            }
        }
        return ButtonEventHandler::mouseEvent(ev);
    }

    void buttonClicked(SubWidget* widget, int button) override
    {
        fIsFocused = true;

        repaint();
        if (fCallback)
            fCallback->focus(this);
    }

    void insertChar(char c) {
        fText = fText.substr(0, cursorPosition) + c + fText.substr(cursorPosition);
        cursorPosition++;
        
        // Schedule text change
        scheduleTextChange(fText);
    }

    void removeChar(void) {
        if (cursorPosition <= 0 || fText.empty()) return;
        
        if (fText[cursorPosition - 1] == ',') {
            std::string toSet = fText.substr(0, cursorPosition - 1);
            bool leftHasPoint = false;
            int i = cursorPosition - 2;
            if (i != -1) {
                while  (i >= 0 && fText[i] != ',') {
                    if (fText[i] == '.') leftHasPoint = true;
                    i--;
                }
            }
            i = cursorPosition;
            while (i < fText.length() && fText[i] != ',') {
                if (((leftHasPoint && fText[i] != '.') || !leftHasPoint) && fText[i] != '-') {
                    toSet += fText[i];
                }
                i++;
            }
            fText = toSet + fText.substr(i, fText.length() - i);
        } else {
            fText = fText.substr(0, cursorPosition - 1) + fText.substr(cursorPosition);
        }
        cursorPosition--;
        
        // Schedule text change
        scheduleTextChange(fText);
    }

    void insertPoint(void) {
        std::string coeff = "";
        int i = cursorPosition - 1;
        while (i >= 0 && fText[i] != ',') {
            coeff = fText[i] + coeff;
            i--;
        }
        i = cursorPosition;
        while (i < fText.length() && fText[i] != ',') {
            coeff += fText[i];
            i++;
        }
        bool found = false;
        for (int i = 0; i < coeff.length(); i++) {
            if (coeff[i] == '.') {
                found = true;
                break;
            }
        }
        if (!found) {
            insertChar('.');
        }
    }

    bool onKeyboard(const KeyboardEvent& ev) override
    {
        if (this->isFocused() && ev.press)
        {
            if (ev.key == 57399 || ev.key == 57397) {
                if (ev.key == 57397 && cursorPosition != 0)
                    cursorPosition--;
                else if (ev.key == 57399 && cursorPosition != fText.length())
                    cursorPosition++;
            } else {
                if (fType == "G" || fType == "M" || fType == " ") {
                    if (ev.key >= '0' && ev.key <= '9')
                        insertChar(ev.key);
                    else if (ev.key == '.')
                        insertPoint();
                    else if (ev.key == '-') {
                        if ((cursorPosition == fText.length() || (fText[cursorPosition] != '-')) && (cursorPosition == 0 || fText[cursorPosition - 1] == ','))
                            insertChar(ev.key);
                    }
                    else if (ev.key == 8 && cursorPosition != 0)
                        removeChar();
                    else if (ev.key == 127 && cursorPosition != fText.length()) {
                        cursorPosition++;
                        removeChar();
                    }
                    else if (ev.key == 13 && fType == "G")
                    {
                        this->setFocused(false);
                        if (fCallback)
                            fCallback->setGain(fText);
                        fText = "";
                        fPoint[fCoeffCount - 1] = false;
                        fMinus[fCoeffCount - 1] = false;
                    }
                    // if (fType == "M" && fCallback) {
                    //     fCallback->setMaxA(fText);
                    // }
                    // if (fType == " " && fCallback) {
                    //     fCallback->generateIfTracking();
                    // }
                } else {
                    if (ev.key >= '0' && ev.key <= '9')
                        insertChar(ev.key);
                    else if (ev.key == ',')
                        insertChar(ev.key);
                    else if (ev.key == 8 && cursorPosition != 0)
                        removeChar();
                    else if (ev.key == 127 && cursorPosition != fText.length()) {
                        cursorPosition++;
                        removeChar();
                    }
                    else if (ev.key == '-') {
                        if ((cursorPosition == fText.length() || (fText[cursorPosition] != '-')) && (cursorPosition == 0 || fText[cursorPosition - 1] == ','))
                            insertChar(ev.key);
                    }
                    else if (ev.key == '.')
                        insertPoint();
                    else if (ev.key == 13)
                        this->setFocused(false);
                    if (fCallback)
                        fCallback->setTF(fType + fText);
                }
            }
            repaint();
        }
        return KeyboardEventHandler::keyboardEvent(ev);
    }

private:
    Callback* fCallback;
    bool fIsFocused;
    bool fPoint[100];
    bool fMinus[100];
    int fCoeffCount;
    std::string fType;
    double textX, textY;
    std::thread fTextChangeThread;
    std::atomic<bool> fTextChangeThreadRunning;
    std::string fPendingTextChange;
    std::mutex fTextChangeMutex;

    void startTextChangeThread() 
    {
        fTextChangeThreadRunning = true;
        fTextChangeThread = std::thread([this]() {
            while (fTextChangeThreadRunning) {
                // Check if there's a pending text change
                if (fTextChanged) {
                    // Check if enough time has passed since the last update
                    auto now = std::chrono::steady_clock::now();
                    auto elapsedMs = std::chrono::duration_cast<std::chrono::milliseconds>(
                        now - fLastTextChangeTime).count();
                        
                    if (elapsedMs >= fTextChangeThrottleMs) {
                        // Get the pending text change
                        std::string textChange;
                        {
                            std::lock_guard<std::mutex> lock(fTextChangeMutex);
                            textChange = fPendingTextChange;
                            fPendingTextChange.clear();
                            fTextChanged = false;
                        }
                        
                        // Send the change to the callback
                        if (fCallback) {
                            if (fType == "M") {
                                fCallback->setMaxA(textChange);
                            } else if (fType == " ") {
                                fCallback->generateIfTracking();
                            } else {
                                fCallback->setTF(fType + textChange);
                            }
                        }
                        
                        fLastTextChangeTime = now;
                    }
                }
                
                // Sleep to avoid consuming too much CPU
                std::this_thread::sleep_for(std::chrono::milliseconds(5));
            }
        });
    }

    void scheduleTextChange(const std::string& text) 
    {
        std::lock_guard<std::mutex> lock(fTextChangeMutex);
        fPendingTextChange = text;
        fTextChanged = true;
    }
};

END_NAMESPACE_DGL
