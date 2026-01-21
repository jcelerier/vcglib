/****************************************************************************
* VCGLib                                                            o o     *
* Visual and Computer Graphics Library                            o     o   *
*                                                                _   O  _   *
* Copyright(C) 2004-2016                                           \/)\/    *
* Visual Computing Lab                                            /\/|      *
* ISTI - Italian National Research Council                           |      *
*                                                                    \      *
* All rights reserved.                                                      *
*                                                                           *
* This program is free software; you can redistribute it and/or modify      *
* it under the terms of the GNU General Public License as published by      *
* the Free Software Foundation; either version 2 of the License, or         *
* (at your option) any later version.                                       *
*                                                                           *
* This program is distributed in the hope that it will be useful,           *
* but WITHOUT ANY WARRANTY; without even the implied warranty of            *
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the             *
* GNU General Public License (http://www.gnu.org/licenses/gpl.txt)          *
* for more details.                                                         *
*                                                                           *
****************************************************************************/
/*! \file trimesh_fast_winding_number.cpp
\ingroup code_sample

\brief Fast winding number evaluation example

This example shows how to use the fast winding number algorithm to evaluate
the generalized winding number at query points on a grid.
*/

#include <vcg/complex/complex.h>
#include <vcg/complex/algorithms/fast_winding_number.h>
#include <vcg/complex/algorithms/create/platonic.h>
#include <wrap/io_trimesh/import_off.h>
#include <wrap/io_trimesh/import_ply.h>
#include <wrap/io_trimesh/export_ply.h>
#include <vcg/complex/algorithms/clean.h>
#include <vcg/complex/algorithms/update/topology.h>
#include <vcg/complex/algorithms/update/normal.h>
#include <vcg/complex/algorithms/update/bounding.h>
#include <vcg/complex/algorithms/update/color.h>
#include <vcg/complex/algorithms/update/quality.h>

class MyVertex;
class MyEdge;
class MyFace;

struct MyUsedTypes : public vcg::UsedTypes<vcg::Use<MyVertex>::AsVertexType,
                                           vcg::Use<MyEdge>::AsEdgeType,
                                           vcg::Use<MyFace>::AsFaceType>
{
};

class MyVertex : public vcg::Vertex<MyUsedTypes, 
                                     vcg::vertex::Coord3f, 
                                     vcg::vertex::Normal3f, 
                                     vcg::vertex::Qualityf,
                                     vcg::vertex::Color4b,
                                     vcg::vertex::BitFlags>
{
};

class MyFace : public vcg::Face<MyUsedTypes, 
                                 vcg::face::FFAdj, 
                                 vcg::face::VertexRef, 
                                 vcg::face::Normal3f, 
                                 vcg::face::BitFlags>
{
};

class MyEdge : public vcg::Edge<MyUsedTypes>
{
};

class MyMesh : public vcg::tri::TriMesh<std::vector<MyVertex>, 
                                         std::vector<MyFace>, 
                                         std::vector<MyEdge>>
{
};

using namespace vcg;

int main(int argc, char **argv)
{
    float accuracyScale = 0.0f;
    int gridRes = 50;
    const char *meshFile = nullptr;
    const char *outputFile = nullptr;
    std::vector<const char*> positionalArgs;

    // Parse command line arguments
    for (int i = 1; i < argc; ++i)
    {
        if (strcmp(argv[i], "--accuracy") == 0)
        {
            if (i + 1 >= argc)
            {
                printf("Error: --accuracy requires a value\n");
                return -1;
            }
            accuracyScale = atof(argv[++i]);
        }
        else if (strcmp(argv[i], "--resolution") == 0)
        {
            if (i + 1 >= argc)
            {
                printf("Error: --resolution requires a value\n");
                return -1;
            }
            gridRes = atoi(argv[++i]);
            if (gridRes < 2)
            {
                printf("Error: grid resolution must be at least 2\n");
                return -1;
            }
        }
        else if (argv[i][0] == '-')
        {
            printf("Error: unknown option %s\n", argv[i]);
            return -1;
        }
        else
        {
            // Positional argument (input file or output file)
            positionalArgs.push_back(argv[i]);
        }
    }

    // Process positional arguments
    if (positionalArgs.empty())
    {
        printf("Usage: %s [options] <input_mesh.off|ply> [output_file.ply]\n", argv[0]);
        printf("\nThis program:\n");
        printf("  - Loads an input mesh\n");
        printf("  - Generates a 2D grid around the mesh\n");
        printf("  - Evaluates the winding number at each grid vertex\n");
        printf("  - Saves the grid with winding numbers as vertex quality\n");
        printf("  - Colors vertices based on winding number (inside/outside)\n");
        printf("\nOptions:\n");
        printf("  --accuracy <value>     : Accuracy threshold (default: 0, exact computation)\n");
        printf("                           Higher values = faster but less accurate\n");
        printf("  --resolution <value>   : Grid resolution NxN (default: 50)\n");
        printf("                           Higher values = more detailed evaluation\n");
        printf("\nArguments:\n");
        printf("  input_mesh             : Input mesh file (OFF or PLY format)\n");
        printf("  output_file            : Optional output file (default: wind_<inputname>.ply)\n");
        return -1;
    }
    
    meshFile = positionalArgs[0];
    
    if (positionalArgs.size() > 1)
    {
        outputFile = positionalArgs[1];
    }
    else
    {
        // Generate default output filename: wind_inputfilename.ply
        std::string inputName(meshFile);
        size_t lastSlash = inputName.find_last_of("/\\");
        if (lastSlash != std::string::npos)
        {
            inputName = inputName.substr(lastSlash + 1);
        }
        // Remove extension
        size_t lastDot = inputName.find_last_of('.');
        if (lastDot != std::string::npos)
        {
            inputName = inputName.substr(0, lastDot);
        }
        static std::string defaultOutput = "wind_" + inputName + ".ply";
        outputFile = defaultOutput.c_str();
    }
    
    if (positionalArgs.size() > 2)
    {
        printf("Error: too many arguments\n");
        return -1;
    }

    printf("Using accuracy scale: %f\n", accuracyScale);
    printf("Using grid resolution: %dx%d\n", gridRes, gridRes);
    printf("Output file: %s\n", outputFile);

    // Load the input mesh
    MyMesh inputMesh;
    int loadmask;
    
    // Try to load as OFF first
    int err = vcg::tri::io::ImporterOFF<MyMesh>::Open(inputMesh, meshFile, loadmask);
    if (err != 0)
    {
        // Try PLY
        err = vcg::tri::io::ImporterPLY<MyMesh>::Open(inputMesh, meshFile, loadmask);
        if (err != 0)
        {
            printf("Error reading file %s\n", meshFile);
            printf("Supported formats: OFF, PLY\n");
            return -1;
        }
    }
    
    printf("Input mesh loaded: %d vertices, %d faces\n", inputMesh.VN(), inputMesh.FN());

    // Prepare input mesh
    vcg::tri::Clean<MyMesh>::RemoveDuplicateVertex(inputMesh);
    vcg::tri::UpdateTopology<MyMesh>::FaceFace(inputMesh);
    vcg::tri::UpdateNormal<MyMesh>::PerFaceNormalized(inputMesh);
    vcg::tri::UpdateBounding<MyMesh>::Box(inputMesh);

    // Initialize fast winding number
    printf("Initializing fast winding number structure...\n");
    tri::FastWindingNumber<MyMesh> fastWN;
    fastWN.init(inputMesh, 2);

    // Create a 2D grid around the mesh
    MyMesh gridMesh;
    Box3f bbox = inputMesh.bbox;
    
    // Expand bbox slightly
    bbox.min *= 1.2;
    bbox.max *= 1.2;
    
    // Grid resolution
    float gridWidth = bbox.DimX();
    float gridHeight = bbox.DimY();
    
    printf("Creating %dx%d grid...\n", gridRes, gridRes);
    
    // Use vcglib built-in Grid function to create a 2D grid mesh
    tri::Grid(gridMesh, gridRes, gridRes, gridWidth, gridHeight);
    
    // Translate grid to match the bounding box position
    for (int i = 0; i < gridMesh.VN(); ++i)
    {
        gridMesh.vert[i].P() += Point3f(bbox.min.X(), bbox.min.Y(), bbox.Center().Z());
    }
    
    printf("Grid created with %d vertices, %d faces\n", gridMesh.VN(), gridMesh.FN());

    // Compute winding numbers for all grid vertices
    printf("Computing winding numbers...\n");
    float minWN = 1e10f, maxWN = -1e10f;
    
    #pragma omp parallel for reduction(min:minWN) reduction(max:maxWN)
    for (int i = 0; i < gridMesh.VN(); ++i)
    {
        Point3f queryPoint = gridMesh.vert[i].P();
        float wn = fastWN.computeSolidAngle(queryPoint, accuracyScale);
        
        // Convert solid angle to winding number (divide by 4*pi)
        wn = wn / (4.0f * M_PI);
        
        gridMesh.vert[i].Q() = wn;
        
        if (wn < minWN) minWN = wn;
        if (wn > maxWN) maxWN = wn;
    }
    
    printf("Winding number range: [%f, %f]\n", minWN, maxWN);

    // Save the grid with winding numbers
    int savemask = vcg::tri::io::Mask::IOM_VERTQUALITY;
    
    if (vcg::tri::io::ExporterPLY<MyMesh>::Save(gridMesh, outputFile, savemask) != 0)
    {
        printf("Error saving output file %s\n", outputFile);
        return -1;
    }
    
    printf("Grid saved to %s\n", outputFile);
    printf("Done!\n");

    return 0;
}
