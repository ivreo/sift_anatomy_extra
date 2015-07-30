/*----------------------------------------------------------------------------

  NLHDR - Non local high dynamic range image generation

  Copyright (c) 2013 cecilia aguerrebere <caguerrebere@gmail.com>
   
  This program is free software: you can redistribute it and/or modify
  it under the terms of the GNU Affero General Public License as
  published by the Free Software Foundation, either version 3 of the
  License, or (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
  GNU Affero General Public License for more details.

  ----------------------------------------------------------------------------*/
  

#include <ImfInputFile.h>
#include <ImfOutputFile.h>
#include <ImfRgbaFile.h>
#include <ImfChannelList.h>
#include <ImfFrameBuffer.h>
#include <half.h>
#include <stdlib.h>
#include <stdio.h>
#include <ImfArray.h>

using namespace Imf;
using namespace Imath;

/* EXR Function Definitions */
struct sRGB
{
  float r, g, b;
};


void writeEXR_float (const char fileName[], const float *zPixels, int width, int height)
{
  Header header (width, height);
  header.channels().insert ("Z", Channel (FLOAT));

  OutputFile file (fileName, header); 
  FrameBuffer frameBuffer;
	
  frameBuffer.insert ("Z", Slice (FLOAT, (char *) zPixels, sizeof (*zPixels) * 1,	sizeof (*zPixels) * width)); 
  file.setFrameBuffer (frameBuffer);
  file.writePixels (height);
}

float *readEXR_float (const char fileName[], int *width_out, int *height_out) {
	
  InputFile file (fileName);
  Array2D<half> rPixels;
  Array2D<half> gPixels;
  Array2D<float> zPixels;
  int width, height;
	
  Box2i dw = file.header().dataWindow();
  width = dw.max.x - dw.min.x + 1;
  height = dw.max.y - dw.min.y + 1;

  rPixels.resizeErase (height, width);
  gPixels.resizeErase (height, width);
  zPixels.resizeErase (height, width);	

  FrameBuffer frameBuffer;
	
  frameBuffer.insert ("R", Slice (HALF, (char *) (&rPixels[0][0] - dw.min.x - dw.min.y * width), sizeof (rPixels[0][0]) * 1, sizeof (rPixels[0][0]) * width, 1, 1, 0.0)); 
  frameBuffer.insert ("G", Slice (HALF, (char *) (&gPixels[0][0] - dw.min.x - dw.min.y * width), sizeof (gPixels[0][0]) * 1, sizeof (gPixels[0][0]) * width, 1, 1, 0.0)); 
  frameBuffer.insert ("Z", Slice (FLOAT,(char *) (&zPixels[0][0] - dw.min.x - dw.min.y * width), sizeof (zPixels[0][0]) * 1, sizeof (zPixels[0][0]) * width, 1, 1, FLT_MAX)); 
	
  file.setFrameBuffer (frameBuffer);
  file.readPixels (dw.min.y, dw.max.y);	

  /* int ydim = height; */
  /* int xdim = width; */
	
  float *zPixels_out = (float*) malloc(height*width*sizeof(float));
  for (int i = 0; i < height; ++i) 
    for (int j = 0; j < width; ++j) 
      zPixels_out[i*width + j] = zPixels[i][j];
		
  *width_out = width;
  *height_out = height;		
	
  return zPixels_out;
}



