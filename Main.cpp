#include <iostream>
#include <map>
#include "Simulation.h"

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

    if (!SDL_GL_SetAttribute(SDL_GL_CONTEXT_MAJOR_VERSION, 1)) return false;
    if (!SDL_GL_SetAttribute(SDL_GL_CONTEXT_MINOR_VERSION, 1)) return false;
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

	while (!quit) {
        handleEvents(quit);
        keys = SDL_GetKeyboardState(NULL);
		if (!handleInput(keys)) quit = true;
		
		Sim->DrawGL();
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

