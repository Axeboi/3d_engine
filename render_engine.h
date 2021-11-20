#pragma once

// TODO
// Matrices needs to be simpler to work with. The matrix init is cumbersome.
// There is a lot of copying of triangles going on each frame. Tranform should take a pointer to tri array instead and re-transform tris back to original position
// Bug: Rasterizer is drawing overlapping pixels. 
// Rasterizer can be rewritten to take advantage of SIMD. This is a good experience to try out.
// Write a loader for obj files. Fun to render something other than cube
// Get camera to work. Needs button controls and camera transform, which should already be in code but commented out. Camera needs looking direction, if it does not already.


#include "olc/olcPixelGameEngine.h"
#include "render_core.h"
#include "render_data_structures.h"
#include "./utils/init_geometry.h"
#include <vector>

#ifndef PI
#define PI    3.14159265358979323846
#endif

double rotation_inc = 0.0;

class RenderEngine : public olc::PixelGameEngine
{
public:
	float f = 0.1;
	float n = 100.0;
	float fov = PI / 2.0;

	TriangleBuffer original_buffer; // Vec4 is 16 bytes. How large should world buffer be?
	TriangleBuffer world_buffer;
	std::vector<std::vector<ScreenBuffer>> screenBuffer;

	float camera_speed{ 0.1f };
	Vec4 camera_location{ 0.0, 0.0, 5.0, 1.0 };
	Vec4 camera_direction{ 0.0, 0.0, -1.0, 1.0 };
	Matrix4 transform;
	Vec4 light {0.0, 0.0, 1.0, 0.0};

	Vec4 move{ 0.0, 0.0, 5.0, 1.0 };

	~RenderEngine()
	{
		triangle_buffer_delete(original_buffer);
  		triangle_buffer_delete(world_buffer);
	}

	RenderEngine()
	{
		sAppName = "3d engine";
	}

public:
	bool OnUserCreate() override
	{
		screenBuffer.resize(ScreenHeight());
		for (int i{ 0 }; i < ScreenHeight(); i++)
		{
			screenBuffer[i].resize(ScreenWidth());
			for (int j{ 0 }; j < ScreenWidth(); j++)
			{
				screenBuffer[i][j] = ScreenBuffer();
			}
		}

		// Initialize buffer and buffer copy
		triangle_buffer_init(original_buffer);

		world_buffer.max_num_tris = original_buffer.max_num_tris;
		world_buffer.size = original_buffer.size;
		world_buffer.tris = (Triangle *) malloc(world_buffer.max_num_tris * sizeof(Triangle));

		//triangle_buffer_init(world_buffer);

		//Vec4 camera {0.0, 0.0, -1.0, 0.0};
		//Vec4 light {0.0, 0.0, -1.0, 0.0};

		//load_scene(original_buffer); 
		return true;
	}

	bool OnUserUpdate(float fElapsedTime) override
	{
		
        //Clear screen
		Clear(olc::BLACK);

		if (GetKey(olc::Key::A).bHeld) camera_location.x -= camera_speed;
		if (GetKey(olc::Key::D).bHeld) camera_location.x += camera_speed;
		if (GetKey(olc::Key::R).bHeld) camera_location.z -= camera_speed;
		if (GetKey(olc::Key::F).bHeld) camera_location.z += camera_speed;
		if (GetKey(olc::Key::W).bHeld) camera_location.y += camera_speed;
		if (GetKey(olc::Key::S).bHeld) camera_location.y -= camera_speed;
		if (GetKey(olc::Key::Z).bPressed) camera_speed = 2.0;
		if (GetKey(olc::Key::X).bPressed) camera_speed = 0.1f;
		//assert_vectors();
		//assert_matrix();

		// Render
		triangle_buffer_copy(&original_buffer, &world_buffer);
		to_identity_matrix(transform);

		//temp_transform_vectors(world_buffer, move);
		get_camera_transform(world_buffer, transform, camera_location, camera_direction);
		apply_light(world_buffer, camera_location);
		get_projection_transform(transform, f, n, fov);
		apply_world_to_projection_transform(world_buffer, transform);
			
		rasterize(world_buffer, screenBuffer, ScreenWidth(), ScreenHeight());
		output();
		reset_screen_buffer(screenBuffer);
		return true;
	}


	void output()
	{
		for (int i = 0; i < ScreenHeight(); i++)
		{
			for (int j = 0; j < ScreenWidth(); j++)
			{
				if (screenBuffer[i][j].r > 0 || screenBuffer[i][j].g > 0 || screenBuffer[i][j].b > 0)
				{
					this->Draw(j, i, olc::Pixel(screenBuffer[i][j].r, screenBuffer[i][j].g, screenBuffer[i][j].b));
				}
			}
		}
	}

	void reset_screen_buffer(std::vector<std::vector<ScreenBuffer>> buf)
	{
		screenBuffer.resize(ScreenHeight());
		for (int i{ 0 }; i < ScreenHeight(); i++)
		{
			screenBuffer[i].resize(ScreenWidth());
			for (int j{ 0 }; j < ScreenWidth(); j++)
			{
				screenBuffer[i][j] = ScreenBuffer();
			}
		}
	}
};
