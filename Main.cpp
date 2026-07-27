#include <iostream>
#include <map>
#include <cstdio>
#include <cstring>
#include "Simulation.h"
#include "VideoRecorder.h"

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
VideoRecorder *Recorder = nullptr;
const bool* keys;
std::map<int, bool> keyPressedMap;
bool SimPaused = false;
double deltaTime = 0.016;

bool init()
{
    if (!SDL_Init(SDL_INIT_VIDEO)) {
        std::cerr << "SDL could not initialize! SDL Error: " << SDL_GetError() << std::endl;
        return false;
    }

    if (!SDL_GL_SetAttribute(SDL_GL_CONTEXT_MAJOR_VERSION, 3)) return false;
    if (!SDL_GL_SetAttribute(SDL_GL_CONTEXT_MINOR_VERSION, 3)) return false;
    if (!SDL_GL_SetAttribute(SDL_GL_CONTEXT_PROFILE_MASK, SDL_GL_CONTEXT_PROFILE_CORE)) return false;
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

    if (!gladLoadGLLoader((GLADloadproc)SDL_GL_GetProcAddress)) {
        std::cerr << "Failed to initialize GLAD" << std::endl;
        return false;
    }

    std::cout << "OpenGL " << glGetString(GL_VERSION) << std::endl;
    std::cout << "GPU: " << glGetString(GL_RENDERER) << std::endl;

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

bool keyPressed(int key, const bool* keys)
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

	double dt = deltaTime;
	if (keys[SDL_SCANCODE_W]) {
		Sim->CamMove(-1.75*dt, 0.0, 0.0);
	}
	if (keys[SDL_SCANCODE_A]) {
		Sim->CamMove( 0.0, 1.75*dt, 0.0);
	}
	if (keys[SDL_SCANCODE_S]) {
		Sim->CamMove( 1.75*dt, 0.0, 0.0);
	}
	if (keys[SDL_SCANCODE_D]) {
		Sim->CamMove( 0.0,-1.75*dt, 0.0);
	}
	if (keys[SDL_SCANCODE_I]) {
		Sim->CamShift(0.0, 420.0*dt, 0.0);
	}
	if (keys[SDL_SCANCODE_K]) {
		Sim->CamShift(0.0, -420.0*dt, 0.0);
	}
	if (keys[SDL_SCANCODE_J]) {
		Sim->CamMove( 0.0, 0.0, -420.0*dt);
	}
	if (keys[SDL_SCANCODE_L]) {
		Sim->CamMove( 0.0, 0.0, 420.0*dt);
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

Uint64 fpsLastTime = 0;
int fpsFrameCount = 0;
double fpsCurrent = 0.0;
Uint64 frameLastTime = 0;

int main(int argc, char *argv[])
{
	std::cout << "---- N-Body Simulation ----" << std::endl;
	std::cout << "Developed by Erick Aldalur" << std::endl;
	std::cout << "2025" << std::endl << std::endl;

	std::string scriptPath = "scripts/default.sim";
	if (argc > 1)
		scriptPath = argv[1];

	bool quit = false;
	SDL_Event e;

	Simulation::ParseDisplaySize(scriptPath, Win.width, Win.height);

	init();

	Sim = new Simulation(scriptPath);
	Sim->ReSizeGL(Win.width, Win.height);

	if (Sim->GetRecordVideo()) {
		std::string videoFile = scriptPath;
		size_t slash = videoFile.find_last_of("/\\");
		std::string scriptDir = "";
		if (slash != std::string::npos) {
			scriptDir = videoFile.substr(0, slash);
			videoFile = videoFile.substr(slash + 1);
		}
		size_t dot = videoFile.rfind('.');
		if (dot != std::string::npos)
			videoFile = videoFile.substr(0, dot);
		videoFile += ".mp4";

		std::string outputDir;
		size_t scriptsPos = scriptDir.find_last_of("/\\");
		if (scriptsPos != std::string::npos)
			outputDir = scriptDir.substr(0, scriptsPos + 1) + "output";
		else
			outputDir = "output";

		SDL_CreateDirectory(outputDir.c_str());
		videoFile = outputDir + "/" + videoFile;
		Recorder = new VideoRecorder(videoFile.c_str(), Win.width, Win.height, 30);
	}

	fpsLastTime = SDL_GetPerformanceCounter();
	frameLastTime = SDL_GetPerformanceCounter();

	while (!quit) {
		Uint64 frameNow = SDL_GetPerformanceCounter();
		deltaTime = (double)(frameNow - frameLastTime) / SDL_GetPerformanceFrequency();
		frameLastTime = frameNow;

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

		if (!Recorder)
			Sim->DrawFPS(fpsCurrent);

		if (Recorder) {
			if (!SimPaused) {
				uint8_t *pixels = new uint8_t[Win.width * Win.height * 3];
				Sim->ReadFramePixels(pixels);
				Recorder->WriteFrame(pixels);
				delete[] pixels;
			} else {
				Sim->BlitToScreen();
			}
		}

		SDL_GL_SwapWindow(Win.gWindow);

		if (!SimPaused)
		{
			if (Sim->GetCamOrbit())
				Sim->CamMove(0.0, Sim->GetCamOrbitTheta(), 0.0);

			Sim->Step();

			if (Sim->GetEndTime() >= 0.0 && Sim->GetTime() >= Sim->GetEndTime())
				quit = true;
		}
    }

	delete Recorder;
	delete Sim;
	close();

	return 0;
}
