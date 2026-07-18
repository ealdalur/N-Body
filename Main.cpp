#include <iostream>
#include <map>
#include <cstdio>
#include <cstring>
#include "Simulation.h"

#ifdef _WIN32
extern "C" {
	__declspec(dllexport) unsigned long NvOptimusEnablement = 1;
	__declspec(dllexport) int AmdPowerXpressRequestHighPerformance = 1;
}
#endif

struct Window
{
	SDL_Window* gWindow = nullptr;
	SDL_GLContext gContext = nullptr;

	int		width,height;
};

Window Win;
Simulation *Sim;
const bool* keys;
std::map<int, bool> keyPressedMap; //for the "keyPressed" function to detect a keypress only once
bool SimPaused = false;

bool init()
{
    if (!SDL_Init(SDL_INIT_VIDEO)) {
        std::cerr << "SDL could not initialize! SDL Error: " << SDL_GetError() << std::endl;
        return false;
    }

    if (!SDL_GL_SetAttribute(SDL_GL_DOUBLEBUFFER, 1)) return false;
    if (!SDL_GL_SetAttribute(SDL_GL_DEPTH_SIZE, 24)) return false;

    Win.gWindow = SDL_CreateWindow(
        "N-Body Simulation",
        Win.width,
        Win.height,
        SDL_WINDOW_OPENGL | SDL_WINDOW_RESIZABLE
    );

    if (Win.gWindow == nullptr) {
        std::cerr << "Window could not be created! SDL Error: " << SDL_GetError() << std::endl;
        return false;
    }

    Win.gContext = SDL_GL_CreateContext(Win.gWindow);
    if (Win.gContext == nullptr) {
        std::cerr << "OpenGL context could not be created! SDL Error: " << SDL_GetError() << std::endl;
        return false;
    }

    SDL_GL_SetSwapInterval(0);

    return true;
}

void handleEvents(bool& quit)
{
    SDL_Event e;
    while (SDL_PollEvent(&e)) {
        if (e.type == SDL_EVENT_QUIT) {
            quit = true;
        }

        if (e.type == SDL_EVENT_WINDOW_RESIZED) {
			Win.width = e.window.data1;
			Win.height = e.window.data2;
            Sim->ReSizeGL(Win.width, Win.height);
        }
    }
}

bool keyPressed(int key, const bool* keys) //this checks if the key is *just* pressed, returns true only once until the key is up again
{
	if (keyPressedMap.find(key) == keyPressedMap.end()) keyPressedMap[key] = false;
	if (keys[key])
	{
		if (keyPressedMap[key] == false)
		{
			keyPressedMap[key] = true;
			return true;
		}
	}
	else keyPressedMap[key] = false;

	return false;
}



bool handleInput(const bool* keys) {

	if (keyPressed(SDL_SCANCODE_G, keys)) {
		//ScreenShot(Win.width,Win.height,"test.bmp");
	}
	if (keyPressed(SDL_SCANCODE_H, keys)) {
		Sim->SaveState();
	}
	if (keyPressed(SDL_SCANCODE_O, keys)) {
		Sim->DrawOctree = !Sim->DrawOctree;
	}
	if (keyPressed(SDL_SCANCODE_M, keys)) {
		Sim->multiThreading = !Sim->multiThreading;
		if (Sim->multiThreading) {
			std::cout << "Number of compute threads: " << Sim->numThreads << std::endl;
		} else {
			std::cout << "Number of compute threads: " << 1 << std::endl;
		}
	}
	if (keyPressed(SDL_SCANCODE_P, keys)) {
		SimPaused = !SimPaused;
		if (SimPaused)
		{
			SDL_SetWindowTitle(Win.gWindow, "N-Body Simulation - Paused");
		}
		else
		{
			SDL_SetWindowTitle(Win.gWindow, "N-Body Simulation");
		}
	}

	if (keys[SDL_SCANCODE_W]) {
		Sim->CamMove(-0.05, 0.0, 0.0);
	}
	if (keys[SDL_SCANCODE_A]) {
		Sim->CamMove( 0.0, 0.05, 0.0);
	}
	if (keys[SDL_SCANCODE_S]) {
		Sim->CamMove( 0.05, 0.0, 0.0);
	}
	if (keys[SDL_SCANCODE_D]) {
		Sim->CamMove( 0.0,-0.05, 0.0);
	}
	if (keys[SDL_SCANCODE_J]) {
		Sim->CamMove( 0.0, 0.0, -12.0);
	}
	if (keys[SDL_SCANCODE_L]) {
		Sim->CamMove( 0.0, 0.0, 12.0);
	}

	if (keys[SDL_SCANCODE_ESCAPE]) {
		return false;
	}
	
	return true;
}

void close()
{
    if (Win.gContext != nullptr) {
        SDL_GL_DestroyContext(Win.gContext);
        Win.gContext = nullptr;
    }

    if (Win.gWindow != nullptr) {
        SDL_DestroyWindow(Win.gWindow);
        Win.gWindow = nullptr;
    }

    SDL_Quit();
}

static const GLubyte font5x7[][7] = {
	{0x00,0x00,0x00,0x00,0x00,0x00,0x00}, // ' ' 32
	{0x04,0x04,0x04,0x04,0x04,0x00,0x04}, // '!' 33
	{0x0A,0x0A,0x00,0x00,0x00,0x00,0x00}, // '"' 34
	{0x0A,0x0A,0x1F,0x0A,0x1F,0x0A,0x0A}, // '#' 35
	{0x04,0x0F,0x14,0x0E,0x05,0x1E,0x04}, // '$' 36
	{0x18,0x19,0x02,0x04,0x08,0x13,0x03}, // '%' 37
	{0x08,0x14,0x14,0x08,0x15,0x12,0x0D}, // '&' 38
	{0x04,0x04,0x00,0x00,0x00,0x00,0x00}, // '\'' 39
	{0x02,0x04,0x08,0x08,0x08,0x04,0x02}, // '(' 40
	{0x08,0x04,0x02,0x02,0x02,0x04,0x08}, // ')' 41
	{0x00,0x04,0x15,0x0E,0x15,0x04,0x00}, // '*' 42
	{0x00,0x04,0x04,0x1F,0x04,0x04,0x00}, // '+' 43
	{0x00,0x00,0x00,0x00,0x00,0x04,0x08}, // ',' 44
	{0x00,0x00,0x00,0x1F,0x00,0x00,0x00}, // '-' 45
	{0x00,0x00,0x00,0x00,0x00,0x00,0x04}, // '.' 46
	{0x00,0x01,0x02,0x04,0x08,0x10,0x00}, // '/' 47
	{0x0E,0x11,0x13,0x15,0x19,0x11,0x0E}, // '0' 48
	{0x04,0x0C,0x04,0x04,0x04,0x04,0x0E}, // '1' 49
	{0x0E,0x11,0x01,0x06,0x08,0x10,0x1F}, // '2' 50
	{0x0E,0x11,0x01,0x06,0x01,0x11,0x0E}, // '3' 51
	{0x02,0x06,0x0A,0x12,0x1F,0x02,0x02}, // '4' 52
	{0x1F,0x10,0x1E,0x01,0x01,0x11,0x0E}, // '5' 53
	{0x06,0x08,0x10,0x1E,0x11,0x11,0x0E}, // '6' 54
	{0x1F,0x01,0x02,0x04,0x08,0x08,0x08}, // '7' 55
	{0x0E,0x11,0x11,0x0E,0x11,0x11,0x0E}, // '8' 56
	{0x0E,0x11,0x11,0x0F,0x01,0x02,0x0C}, // '9' 57
	{0x00,0x00,0x04,0x00,0x04,0x00,0x00}, // ':' 58
	{0x00,0x00,0x04,0x00,0x04,0x04,0x08}, // ';' 59
	{0x02,0x04,0x08,0x10,0x08,0x04,0x02}, // '<' 60
	{0x00,0x00,0x1F,0x00,0x1F,0x00,0x00}, // '=' 61
	{0x08,0x04,0x02,0x01,0x02,0x04,0x08}, // '>' 62
	{0x0E,0x11,0x01,0x02,0x04,0x00,0x04}, // '?' 63
	{0x0E,0x11,0x17,0x15,0x17,0x10,0x0E}, // '@' 64
	{0x0E,0x11,0x11,0x1F,0x11,0x11,0x11}, // 'A' 65
	{0x1E,0x11,0x11,0x1E,0x11,0x11,0x1E}, // 'B' 66
	{0x0E,0x11,0x10,0x10,0x10,0x11,0x0E}, // 'C' 67
	{0x1E,0x11,0x11,0x11,0x11,0x11,0x1E}, // 'D' 68
	{0x1F,0x10,0x10,0x1E,0x10,0x10,0x1F}, // 'E' 69
	{0x1F,0x10,0x10,0x1E,0x10,0x10,0x10}, // 'F' 70
	{0x0E,0x11,0x10,0x17,0x11,0x11,0x0F}, // 'G' 71
	{0x11,0x11,0x11,0x1F,0x11,0x11,0x11}, // 'H' 72
	{0x0E,0x04,0x04,0x04,0x04,0x04,0x0E}, // 'I' 73
	{0x07,0x02,0x02,0x02,0x02,0x12,0x0C}, // 'J' 74
	{0x11,0x12,0x14,0x18,0x14,0x12,0x11}, // 'K' 75
	{0x10,0x10,0x10,0x10,0x10,0x10,0x1F}, // 'L' 76
	{0x11,0x1B,0x15,0x15,0x11,0x11,0x11}, // 'M' 77
	{0x11,0x19,0x15,0x13,0x11,0x11,0x11}, // 'N' 78
	{0x0E,0x11,0x11,0x11,0x11,0x11,0x0E}, // 'O' 79
	{0x1E,0x11,0x11,0x1E,0x10,0x10,0x10}, // 'P' 80
	{0x0E,0x11,0x11,0x11,0x15,0x12,0x0D}, // 'Q' 81
	{0x1E,0x11,0x11,0x1E,0x14,0x12,0x11}, // 'R' 82
	{0x0E,0x11,0x10,0x0E,0x01,0x11,0x0E}, // 'S' 83
	{0x1F,0x04,0x04,0x04,0x04,0x04,0x04}, // 'T' 84
	{0x11,0x11,0x11,0x11,0x11,0x11,0x0E}, // 'U' 85
	{0x11,0x11,0x11,0x11,0x11,0x0A,0x04}, // 'V' 86
	{0x11,0x11,0x11,0x15,0x15,0x1B,0x11}, // 'W' 87
	{0x11,0x11,0x0A,0x04,0x0A,0x11,0x11}, // 'X' 88
	{0x11,0x11,0x0A,0x04,0x04,0x04,0x04}, // 'Y' 89
	{0x1F,0x01,0x02,0x04,0x08,0x10,0x1F}, // 'Z' 90
};

void drawChar(char c, float x, float y, int scale)
{
	if (c < 32 || c > 90) return;
	const GLubyte *glyph = font5x7[c - 32];

	for (int row = 0; row < 7; row++) {
		for (int col = 0; col < 5; col++) {
			if (glyph[6 - row] & (1 << (4 - col))) {
				glVertex2f(x + col * scale, y + row * scale);
				glVertex2f(x + col * scale + scale, y + row * scale);
				glVertex2f(x + col * scale + scale, y + row * scale + scale);
				glVertex2f(x + col * scale, y + row * scale + scale);
			}
		}
	}
}

void drawString(const char *str, float x, float y, int scale)
{
	glMatrixMode(GL_PROJECTION);
	glPushMatrix();
	glLoadIdentity();
	glOrtho(0, Win.width, 0, Win.height, -1, 1);
	glMatrixMode(GL_MODELVIEW);
	glPushMatrix();
	glLoadIdentity();

	glColor4f(1.0f, 1.0f, 1.0f, 0.8f);
	glBegin(GL_QUADS);
	for (const char *p = str; *p; p++) {
		drawChar(*p, x, y, scale);
		x += 6 * scale;
	}
	glEnd();

	glPopMatrix();
	glMatrixMode(GL_PROJECTION);
	glPopMatrix();
	glMatrixMode(GL_MODELVIEW);
}

Uint64 fpsLastTime = 0;
int fpsFrameCount = 0;
double fpsCurrent = 0.0;

int main()
{
	std::cout << "---- N-Body Simulation ----" << std::endl;
	std::cout << "Developed by Erick Aldalur" << std::endl;
	std::cout << "2025" << std::endl << std::endl; 

	bool quit = false;
	SDL_Event e;

	Win.width = 2*1920/3;
	Win.height = 2*1080/3;

	init();
	
	Sim = new Simulation();
	Sim->ReSizeGL(Win.width, Win.height);

	fpsLastTime = SDL_GetPerformanceCounter();

	while (!quit) {
        handleEvents(quit);
        keys = SDL_GetKeyboardState(NULL);
		if (!handleInput(keys)) quit = true;

		Sim->DrawGL();

		fpsFrameCount++;
		Uint64 now = SDL_GetPerformanceCounter();
		double elapsed = (double)(now - fpsLastTime) / SDL_GetPerformanceFrequency();
		if (elapsed >= 0.5) {
			fpsCurrent = fpsFrameCount / elapsed;
			fpsFrameCount = 0;
			fpsLastTime = now;
		}

		char fpsText[16];
		snprintf(fpsText, sizeof(fpsText), "FPS:%d", (int)(fpsCurrent + 0.5));
		int textWidth = (int)strlen(fpsText) * 6;
		drawString(fpsText, Win.width - textWidth - 4, Win.height - 12, 1);

		SDL_GL_SwapWindow(Win.gWindow);

		if (!SimPaused)
		{
			Sim->Step();
		}
    }

	delete Sim;
	close();

	return 0;
}

