#pragma once

#include <string>
#include <fstream>
#include <sstream>
#include <iostream>
#include <assert.h>
#include <vector>
#include "../render_data_structures.h"

const std::string file_name = "C:\\Users\\fredr\\source\\repos\\graphics\\3dengine\\assets\\mario_intro_stage.obj";
//const std::string file_name = "C:\\Users\\fredr\\source\\repos\\graphics\\3dengine\\assets\\test.obj";
//const std::string file_name = "../assets/mario_intro_stage.obj";

void parse_obj_file(TriangleBuffer &);

void
triangle_buffer_init(TriangleBuffer &buffer)
{
    //std::string file_name("../assets/mario_intro_stage.obj");
    std::string line;
    std::fstream obj_file(file_name);
    if(!obj_file)
    {
        std::cout << "File load error: could not load " << file_name << std::endl;
        exit(-1);
    }

    uint32_t counter = 0;
    while(std::getline(obj_file, line))
    {
        if (line.substr(0,2) == "v ")
            counter++;
    }

    //assert(counter % 3 == 0);

    buffer.max_num_tris = 100000;
    buffer.size = 0;
    buffer.tris = (Triangle *) malloc(buffer.max_num_tris * sizeof(Triangle));
    if(buffer.tris == NULL)
    {
        std::cout << "Allocation error: malloc failed. Exiting...";
        exit(-1);
    }

    parse_obj_file(buffer);
}

// void
// parse_obj_file(TriangleBuffer &buffer)
// {
//     //std::string file_name("../assets/mario_intro_stage.obj");
//     std::string line;
//     std::fstream obj_file(file_name);

//     uint32_t buffer_index = 0;
//     while(std::getline(obj_file, line))
//     {
//         if(line.substr(0, 2) != "v ")
//             continue;
        
//         Triangle triangle;
//         for (size_t i = 0; i < 3; i++)
//         {
//             std::istringstream vector(line.substr(2));
//             vector >> triangle.tri[i].x;
//             vector >> triangle.tri[i].y;
//             vector >> triangle.tri[i].z;
//         }

//         buffer.tris[buffer_index] = triangle;
//         buffer.size++;
        
//         buffer_index++;
//     }
// }

void
parse_obj_file(TriangleBuffer &buffer)
{
    //std::string file_name("../assets/mario_intro_stage.obj");
    std::string line;
    std::fstream obj_file(file_name);

    uint32_t buffer_index = 0;
    std::vector<Vec4> vertex_arr;
    std::vector<Vec4> normal_arr;

    std::string face_delimiter = "/";
    while(std::getline(obj_file, line))
    {
        if(line.substr(0, 2) == "v ")
        { 
            Vec4 point_to_array;
            std::istringstream point_from_file(line.substr(2));
            point_from_file >> point_to_array.x;
            point_from_file >> point_to_array.y;
            point_from_file >> point_to_array.z;
            point_to_array.w = 1.0;
            vertex_arr.push_back(point_to_array);
        }

        if(line.substr(0,2) == "vt")
        {

        }

        if(line.substr(0,2) == "vn")
        {
            Vec4 normal_to_array;
            std::istringstream normal_from_file(line.substr(2));
            normal_from_file >> normal_to_array.x;
            normal_from_file >> normal_to_array.y;
            normal_from_file >> normal_to_array.z;
            normal_to_array.w = 1.0;

            // normal_to_array.x = -normal_to_array.x;
            // normal_to_array.y = -normal_to_array.y;
            // normal_to_array.z = -normal_to_array.z;

            normal_arr.push_back(normal_to_array);
        }

        if(line.substr(0, 2) == "f ")
        {
            Triangle triangle;
            std::istringstream faces_from_file(line.substr(2));
            int i = 2;
            std::string part;
            while(faces_from_file >> part)
            {
                size_t start = 0;
                size_t end = 0;
                std::string token;

                end = part.find(face_delimiter, start);
                token = part.substr(0, end);

                triangle.tri[i] = vertex_arr[stoi(token) - 1];
                start = end + 1;

                end = part.find(face_delimiter, start);
                token = part.substr(start, end);
                start = end + 1;

                token = part.substr(start);

                triangle.normal = normal_arr[stoi(token) - 1];
            
                i--;
            }

            buffer.tris[buffer_index] = triangle;
            buffer.size++;
            buffer_index++;
        }        
    }
}
