#include <windows.h>
#include <stdio.h>
#include "Simulation.h"

struct Window
{
	HDC			hDC;
	HGLRC		hRC;
	HWND		hWnd;
	HINSTANCE	hInstance;

	int		width,height;

	bool	keys[256];
	bool	active;
};

Window Win;
Simulation *Sim;
bool SimPaused = false;
bool gp = false;
bool hp = false;
bool op = false;
bool pp = false;

void ScreenShot(int windowWidth, int windowHeight, char* filename)
{
	byte* Buff = new byte[windowWidth*windowHeight*3];
	byte temp;
	if (!Buff)
		return;

	glReadBuffer(GL_BACK);
	glReadPixels(0, 0, windowWidth, windowHeight, GL_RGB, GL_UNSIGNED_BYTE, Buff);

	int x_off,y_off;
	for (int i=0; i<windowWidth; i++)
	{
		x_off = 3*i;
		for (int j =0; j<windowHeight; j++)
		{
			y_off = 3*windowWidth*j;
			
			temp = Buff[y_off+x_off+0];
			Buff[y_off+x_off+0] = Buff[y_off+x_off+2];
			Buff[y_off+x_off+2] = temp;
		}
	}

	FILE *Out;
	fopen_s(&Out, filename, "wb");

	if (!Out)
		return;
	BITMAPFILEHEADER bitmapFileHeader;
	BITMAPINFOHEADER bitmapInfoHeader;

	bitmapFileHeader.bfType = 0x4D42;
	bitmapFileHeader.bfSize = windowWidth*windowHeight * 3;
	bitmapFileHeader.bfReserved1 = 0;
	bitmapFileHeader.bfReserved2 = 0;
	bitmapFileHeader.bfOffBits =
		sizeof(BITMAPFILEHEADER) + sizeof(BITMAPINFOHEADER);

	bitmapInfoHeader.biSize = sizeof(BITMAPINFOHEADER);
	bitmapInfoHeader.biWidth = windowWidth;
	bitmapInfoHeader.biHeight = windowHeight;
	bitmapInfoHeader.biPlanes = 1;
	bitmapInfoHeader.biBitCount = 24;
	bitmapInfoHeader.biCompression = BI_RGB;
	bitmapInfoHeader.biSizeImage = 0;
	bitmapInfoHeader.biXPelsPerMeter = 0; // ?
	bitmapInfoHeader.biYPelsPerMeter = 0; // ?
	bitmapInfoHeader.biClrUsed = 0;
	bitmapInfoHeader.biClrImportant = 0;

	fwrite(&bitmapFileHeader, sizeof(BITMAPFILEHEADER), 1, Out);
	fwrite(&bitmapInfoHeader, sizeof(BITMAPINFOHEADER), 1, Out);
	fwrite(Buff, windowWidth*windowHeight * 3, 1, Out);
	fclose(Out);

	delete[] Buff;
}

void ProcessControls()
{
	if (Win.keys['G'] && (!gp))
	{
		gp = true;
		ScreenShot(Win.width,Win.height,"test.bmp");
	}
	if (gp && (!Win.keys['G']))
	{
		gp = false;
	}
	
	if (Win.keys['H'] && (!hp))
	{
		hp = true;
		Sim->SaveState();
	}
	if (hp && (!Win.keys['H']))
	{
		hp = false;
	}
	
	if (Win.keys['O'] && (!op))
	{
		op = true;
		Sim->DrawOctree = !Sim->DrawOctree;
	}
	if (op && (!Win.keys['O']))
	{
		op = false;
	}
	if (Win.keys['P'] && (!pp))
	{
		pp = true;
		SimPaused = !SimPaused;
		if (SimPaused)
		{
			SetWindowText(Win.hWnd, "N-Body Simulation - Paused");
		}
		else
		{
			SetWindowText(Win.hWnd, "N-Body Simulation");
		}
	}
	if (pp && (!Win.keys['P']))
	{
		pp = false;
	}

	if (Win.keys['W']) 
	{
		Sim->CamMove(-0.05, 0.0, 0.0);
	}
	if (Win.keys['A']) 
	{
		Sim->CamMove( 0.0, 0.05, 0.0);
	}
	if (Win.keys['S']) 
	{
		Sim->CamMove( 0.05, 0.0, 0.0);
	}
	if (Win.keys['D']) 
	{
		Sim->CamMove( 0.0,-0.05, 0.0);
	}
	if (Win.keys['J']) 
	{
		Sim->CamMove( 0.0, 0.0, -12.0);
	}
	if (Win.keys['L']) 
	{
		Sim->CamMove( 0.0, 0.0, 12.0);
	}
	
}

LRESULT CALLBACK WndProc(HWND hWnd, UINT uMsg, WPARAM wParam, LPARAM lParam)
{
	switch (uMsg)
	{
		case WM_ACTIVATE:
		{
			// LoWord Can Be WA_INACTIVE, WA_ACTIVE, WA_CLICKACTIVE,
			// The High-Order Word Specifies The Minimized State Of The Window Being Activated Or Deactivated.
			// A NonZero Value Indicates The Window Is Minimized.
			if ((LOWORD(wParam) != WA_INACTIVE) && !((BOOL)HIWORD(wParam)))
				Win.active=TRUE;
			else
				Win.active=FALSE;

			return 0;
		}

		case WM_SYSCOMMAND:
		{
			switch (wParam)	
			{
				case SC_SCREENSAVE:
				case SC_MONITORPOWER:
				return 0;
			}
			break;
		}

		case WM_CLOSE:
		{
			PostQuitMessage(0);
			return 0;
		}

		case WM_KEYDOWN:
		{
			Win.keys[wParam] = TRUE;
			return 0;
		}

		case WM_KEYUP:
		{
			Win.keys[wParam] = FALSE;
			return 0;
		}

		case WM_SIZE:
		{
			Win.width = LOWORD(lParam);
			Win.height = HIWORD(lParam);
			Sim->ReSizeGL(Win.width,Win.height);
			return 0;
		}
	}

	return DefWindowProc(hWnd,uMsg,wParam,lParam);
}

void CreateGLWindow(LPCSTR Name, int width, int height, BYTE bits)
{
	GLuint		PixelFormat;
	WNDCLASS	wc;
	DWORD		dwExStyle;
	DWORD		dwStyle;
	RECT		WindowRect;
	WindowRect.left=(long)0;
	WindowRect.right=(long)width;
	WindowRect.top=(long)0;
	WindowRect.bottom=(long)height;
	
	Win.width = width;
	Win.height = height;

	Win.hInstance		= GetModuleHandle(NULL);
	wc.style			= CS_HREDRAW | CS_VREDRAW | CS_OWNDC;
	wc.lpfnWndProc		= (WNDPROC)WndProc;
	wc.cbClsExtra		= 0;
	wc.cbWndExtra		= 0;
	wc.hInstance		= Win.hInstance;
	wc.hIcon			= LoadIcon(NULL, IDI_WINLOGO);
	wc.hCursor			= LoadCursor(NULL, IDC_ARROW);
	wc.hbrBackground	= NULL;
	wc.lpszMenuName		= NULL;
	wc.lpszClassName	= "GLWindow";

	RegisterClass(&wc);

	dwExStyle=WS_EX_APPWINDOW | WS_EX_WINDOWEDGE;
	dwStyle=WS_OVERLAPPEDWINDOW;

	AdjustWindowRectEx(&WindowRect, dwStyle, FALSE, dwExStyle);

	// Create The Window
	Win.hWnd=CreateWindowEx(dwExStyle,							// Extended Style For The Window
							"GLWindow",							// Class Name
							Name,								// Window Title
							dwStyle |							// Defined Window Style
							WS_CLIPSIBLINGS |					// Required Window Style
							WS_CLIPCHILDREN,					// Required Window Style
							100, 100,							// Window Position
							WindowRect.right-WindowRect.left,	// Calculate Window Width
							WindowRect.bottom-WindowRect.top,	// Calculate Window Height
							NULL,								// No Parent Window
							NULL,								// No Menu
							Win.hInstance,						// Instance
							NULL);								// Dont Pass Anything To WM_CREATE

	static	PIXELFORMATDESCRIPTOR pfd=
	{
		sizeof(PIXELFORMATDESCRIPTOR),				// Size Of This Pixel Format Descriptor
		1,											// Version Number
		PFD_DRAW_TO_WINDOW |						// Format Must Support Window
		PFD_SUPPORT_OPENGL |						// Format Must Support OpenGL
		PFD_DOUBLEBUFFER,							// Must Support Double Buffering
		PFD_TYPE_RGBA,								// Request An RGBA Format
		bits,										// Select Our Color Depth
		0, 0, 0, 0, 0, 0,							// Color Bits Ignored
		0,											// No Alpha Buffer
		0,											// Shift Bit Ignored
		0,											// No Accumulation Buffer
		0, 0, 0, 0,									// Accumulation Bits Ignored
		16,											// 16Bit Z-Buffer (Depth Buffer)  
		0,											// No Stencil Buffer
		0,											// No Auxiliary Buffer
		PFD_MAIN_PLANE,								// Main Drawing Layer
		0,											// Reserved
		0, 0, 0										// Layer Masks Ignored
	};
	
	Win.hDC=GetDC(Win.hWnd);
	
	PixelFormat=ChoosePixelFormat(Win.hDC,&pfd);
	SetPixelFormat(Win.hDC,PixelFormat,&pfd);
	
	Win.hRC=wglCreateContext(Win.hDC);
	wglMakeCurrent(Win.hDC,Win.hRC);

	ShowWindow(Win.hWnd,SW_SHOW);
	SetForegroundWindow(Win.hWnd);
	SetFocus(Win.hWnd);
}

int main()
{
	MSG		msg;
	BOOL	done=FALSE;
	DWORD t = GetTickCount();

	Win.width = 2*1920/3;
	Win.height = 2*1080/3;

	CreateGLWindow("N-Body Simulation",Win.width,Win.height,16);
	
	Sim = new Simulation();
	Sim->ReSizeGL(Win.width, Win.height);

	while(!done)
	{
		if (PeekMessage(&msg,NULL,0,0,PM_REMOVE))
		{
			if (msg.message==WM_QUIT)
			{
				done=TRUE;
			}
			else
			{
				TranslateMessage(&msg);
				DispatchMessage(&msg);
			}
		}
		else
		{
			if (Win.keys[VK_ESCAPE])
			{
				done=TRUE;
			}
			else
			{
				if ((GetTickCount()-t) > 10)
				{
					ProcessControls();
					Sim->DrawGL();
					
					if (!SimPaused)
					{
						Sim->Step();
					}

					glFlush();
					SwapBuffers(Win.hDC);

					t = GetTickCount();
				}
			}
		}
	}

	delete Sim;

	wglMakeCurrent(NULL,NULL);
	wglDeleteContext(Win.hRC);
	ReleaseDC(Win.hWnd,Win.hDC);
	DestroyWindow(Win.hWnd);
	UnregisterClass("GLWindow",Win.hInstance);

	return 0;
}

