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
  

#ifndef IO_EXR_H_
#define IO_EXR_H_


void writeEXR_float(const char fileName[], const float *zPixels, int width, int height);
float *readEXR_float (const char fileName[], int *width_out, int *height_out);

#endif
