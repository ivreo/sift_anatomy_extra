# Copyright 2014 Ives Rey-Otero <ivesreyotero@gmail.com>

# compilers configuration
#CC = gcc
CC = clang
OFLAGS = -g -O3
LIBS = -lpng -L/usr/local/lib -lm
#-lquadmath   -ltcmalloc

CFLAGS = -Wall -Wno-write-strings -pedantic -std=c99 -D_POSIX_C_SOURCE=200809L  #   -Werror

# Source files with executables.
SRC_ALGO = sift_cli

SRC_MATCH = match_cli

# TEMP normalized_patch
SRC_DEMO = #demo_extract_patch

#SRCa = lib_sift.c \
#	   lib_sift_anatomy.c \
#	   lib_scalespace.c \
#	   lib_description.c \
#       lib_discrete.c \
#	   lib_keypoint.c \
#	   lib_util.c\
#	   lib_ellipse_anatomy.c

SRCa = lib_sift_anatomy.c \
	   lib_scalespace.c \
	   lib_description.c \
       lib_discrete.c \
	   lib_keypoint.c \
	   lib_bsplines.c \
	   lib_util.c

SRCb = lib_io_scalespace.c

SRCc = lib_matching.c

SRCDIR = src
OBJDIR = src
BINDIR = bin

OBJa = $(addprefix $(OBJDIR)/,$(SRCa:.c=.o))
OBJb = $(addprefix $(OBJDIR)/,$(SRCb:.c=.o))
OBJc = $(addprefix $(OBJDIR)/,$(SRCc:.c=.o))

OBJ = $(OBJa) $(OBJb) $(OBJc)

BIN = $(addprefix $(BINDIR)/,$(SRC_ALGO))
BINMATCH = $(addprefix $(BINDIR)/,$(SRC_MATCH))
BINDEMO = $(addprefix $(BINDIR)/,$(SRC_DEMO))

sift= $(BIN)
match= $(BINMATCH)
demo= $(BINDEMO)
#default: $(OBJDIR) $(BINDIR) $(sift) $(match) $(demo)
#extra= $(BINDIR) $(BINDIR)/match_cli_ellipse $(BINDIR)/normalized_patch $(BINDIR)/find_orientation_and_describe_keys $(BINDIR)/extra_sift_cli $(BINDIR)/gradual
extra= $(BINDIR) $(BINDIR)/gradual    $(BINDIR)/write_dct_scalespace $(BINDIR)/crop_from_scalespace  $(BINDIR)/extra_sift_cli   $(BINDIR)/grab_ss_neighborhood $(BINDIR)/check_extremum_in_cube $(BINDIR)/descriptors_from_scalespace  $(BINDIR)/check_extrema_in_scalespace $(BINDIR)/check_extrema_from_scalespace $(BINDIR)/detection_from_scalespace  #$(BINDIR)/extra_sift_cli_extrema_patches
default: $(OBJDIR) $(BINDIR) $(sift) $(match) $(demo) $(extra)

#---------------------------------------------------------------
#  SIFT CLI
#

$(BIN) : $(BINDIR)/% : $(SRCDIR)/%.c $(OBJDIR)/lib_sift_anatomy.o  $(OBJDIR)/lib_keypoint.o $(OBJDIR)/lib_scalespace.o $(OBJDIR)/lib_description.o   $(OBJDIR)/lib_discrete.o  $(OBJDIR)/lib_io_scalespace.o  $(OBJDIR)/lib_util.o   $(OBJDIR)/io_png.o  $(OBJDIR)/lib_bsplines.o  $(OBJDIR)/iio.o
	$(CC) $(CFLAGS) -o $@ $^ $(LIBS) -ltiff -ljpeg
	#$(OBJDIR)/lib_sift.o


$(OBJDIR):
	    -mkdir -p $(OBJDIR)

$(BINDIR):
	    -mkdir -p $(BINDIR)

#---------------------------------------------------------------
#  LIB_SIFT
#
#$(OBJDIR)/lib_sift.o : $(SRCDIR)/lib_sift.c $(OBJDIR)/lib_sift_anatomy.o $(OBJDIR)/lib_keypoint.o $(OBJDIR)/lib_util.o
#	$(CC) $(CFLAGS) $(OFLAGS) -c -o $@ $<

$(OBJDIR)/lib_scalespace.o : $(SRCDIR)/lib_scalespace.c $(OBJDIR)/lib_discrete.o  $(OBJDIR)/lib_util.o 
	$(CC) $(CFLAGS) $(OFLAGS) -c -o $@ $<

$(OBJDIR)/lib_discrete.o : $(SRCDIR)/lib_discrete.c $(OBJDIR)/lib_util.o 
	$(CC) $(CFLAGS) $(OFLAGS) -c -o $@ $<

$(OBJDIR)/lib_description.o : $(SRCDIR)/lib_description.c $(OBJDIR)/lib_discrete.o $(OBJDIR)/lib_keypoint.o $(OBJDIR)/lib_util.o 
	$(CC) $(CFLAGS) $(OFLAGS) -c -o $@ $<

$(OBJDIR)/lib_keypoint.o : $(SRCDIR)/lib_keypoint.c $(OBJDIR)/lib_util.o 
	$(CC) $(CFLAGS) $(OFLAGS) -c -o $@ $<

$(OBJDIR)/lib_sift_anatomy.o : $(SRCDIR)/lib_sift_anatomy.c $(OBJDIR)/lib_keypoint.o $(OBJDIR)/lib_discrete.o $(OBJDIR)/lib_scalespace.o $(OBJDIR)/lib_util.o 
	$(CC) $(CFLAGS) $(OFLAGS) -c -o $@ $<

$(OBJDIR)/lib_util.o : $(SRCDIR)/lib_util.c
	$(CC) $(CFLAGS) $(OFLAGS) -c -o $@ $<

#--------------------------------------------------------------
#   IN (image) and OUT (scalespace)
#
$(OBJDIR)/io_png.o : $(SRCDIR)/io_png.c
	$(CC) $(CFLAGS) $(OFLAGS) -c -o $@ $<

$(OBJDIR)/lib_io_scalespace.o : $(SRCDIR)/lib_io_scalespace.c $(OBJDIR)/io_png.o  $(OBJDIR)/lib_scalespace.o $(OBJDIR)/lib_util.o
	$(CC) $(CFLAGS) $(OFLAGS) -c -o $@ $<


#-------------------------------------------------------------
#   Matching algorithm
#
$(OBJDIR)/lib_matching.o : $(SRCDIR)/lib_matching.c $(OBJDIR)/lib_keypoint.o $(OBJDIR)/lib_util.o
	$(CC) $(CFLAGS) $(OFLAGS) -c -o $@ $<

$(BINMATCH) : $(SRCDIR)/match_cli.c $(OBJDIR)/lib_keypoint.o $(OBJDIR)/lib_matching.o   $(OBJDIR)/lib_util.o
	$(CC) $(CFLAGS) -o $@ $^ $(LIBS) -lm

#-------------------------------------------------------------
#  Tools used in the demo 
#
$(BINDEMO) : $(BINDIR)/% :	 $(SRCDIR)/demo_extract_patch.c  $(OBJDIR)/lib_discrete.o $(OBJDIR)/io_png.o $(OBJDIR)/lib_util.o
	$(CC) $(CFLAGS) -o $@ $^ $(LIBS)



#---------------------- EXTRA EXTRA EXTRA EXTRA ------------------------------
#---------------------- EXTRA EXTRA EXTRA EXTRA ------------------------------
#---------------------- EXTRA EXTRA EXTRA EXTRA ------------------------------





#------------- ELLIPSE---------------------------------------------
#$(OBJDIR)/lib_ellipse_anatomy.o : $(SRCDIR)/lib_ellipse_anatomy.c
#	$(CC) $(CFLAGS) $(OFLAGS) -c -o $@ $^

# ELLIPSE MATCHING
#$(BINDIR)/match_cli_ellipse : $(SRCDIR)/match_cli_ellipse.c $(OBJ) $(OBJDIR)/io_png.o $(OBJDIR)/lib_ellipse_anatomy.o
#	$(CC) $(CFLAGS) $(OFLAGS) -o $@ $^ $(LIBS)

# AFFINE NORMALIZED PATCH + DESCRIPTION (ORI + FEAT)
$(BINDIR)/normalized_patch : $(SRCDIR)/normalized_patch.c $(OBJ) $(OBJDIR)/io_png.o $(OBJDIR)/lib_ellipse_anatomy.o
	$(CC) $(CFLAGS) $(OFLAGS) -o $@ $^ $(LIBS)

# DESCRIPTION (ORI + FEAT)
$(BINDIR)/find_orientation_and_describe_keys : $(SRCDIR)/find_orientation_and_describe_keys.c $(OBJ) $(OBJDIR)/io_png.o

#------------- EXR + IIO ---------------------------------------------
CFLAGS2 = -Wno-write-strings -std=c99 -D_POSIX_C_SOURCE=200809L # -pedantic -Werror
$(OBJDIR)/iio.o :  $(SRCDIR)/iio.c
	$(CC) $(CFLAGS2) $(OFLAGS) -c -o $@ $^

#CPP = g++
CPP = clang++

CPPFLAGS = -Wall -Wno-write-strings -D_POSIX_C_SOURCE=200809L  # -pedantic -Werror
$(OBJDIR)/io_exr.o : $(SRCDIR)/io_exr.cpp
	$(CPP) $(CPPFLAGS) $(OFLAGS) -c -o $@ $^ -I/usr/include/OpenEXR

$(OBJDIR)/lib_fourier.o : $(SRCDIR)/lib_fourier.c
	$(CC) $(CFLAGS) $(OFLAGS) -c -o $@ $^

$(OBJDIR)/lib_bsplines.o : $(SRCDIR)/lib_bsplines.c
	$(CC) $(CFLAGS) $(OFLAGS) -c -o $@ $^

$(OBJDIR)/lib_discrete_extra.o : $(SRCDIR)/lib_discrete_extra.c
	$(CC) $(CFLAGS) $(OFLAGS) -c -o $@ $^

$(OBJDIR)/lib_dense_anatomy.o : $(SRCDIR)/lib_dense_anatomy.c
	$(CC) $(CFLAGS) $(OFLAGS) -c -o $@ $^

$(BINDIR)/extra_sift_cli : $(SRCDIR)/extra_sift_cli.cpp $(OBJ) $(OBJDIR)/iio.o $(OBJDIR)/io_png.o  $(OBJDIR)/io_exr.o  $(OBJDIR)/lib_fourier.o $(OBJDIR)/lib_dense_anatomy.o $(OBJDIR)/lib_discrete_extra.o
	$(CPP) $(CPPFLAGS) $(OFLAGS) -o $@ $^ $(LIBS) -ltiff -ljpeg -lIlmImf -lHalf -lfftw3 -lfftw3f #-lfftw3q

$(BINDIR)/gradual : $(SRCDIR)/gradual.cpp $(OBJ) $(OBJDIR)/iio.o $(OBJDIR)/io_png.o  $(OBJDIR)/io_exr.o  $(OBJDIR)/lib_fourier.o $(OBJDIR)/lib_dense_anatomy.o $(OBJDIR)/lib_discrete_extra.o
	$(CPP) $(CPPFLAGS) $(OFLAGS) -o $@ $^ $(LIBS) -ltiff -ljpeg -lIlmImf -lHalf -lfftw3 -lfftw3f #-lfftw3q

$(BINDIR)/write_dct_scalespace : $(SRCDIR)/write_dct_scalespace.cpp $(OBJ) $(OBJDIR)/iio.o $(OBJDIR)/io_png.o  $(OBJDIR)/io_exr.o  $(OBJDIR)/lib_fourier.o $(OBJDIR)/lib_dense_anatomy.o $(OBJDIR)/lib_discrete_extra.o
	$(CPP) $(CPPFLAGS) $(OFLAGS) -o $@ $^ $(LIBS) -ltiff -ljpeg -lIlmImf -lHalf -lfftw3 -lfftw3f #-lfftw3q

$(BINDIR)/crop_from_scalespace : $(SRCDIR)/crop_from_scalespace.cpp $(OBJ) $(OBJDIR)/iio.o $(OBJDIR)/io_png.o  $(OBJDIR)/io_exr.o  $(OBJDIR)/lib_fourier.o $(OBJDIR)/lib_dense_anatomy.o $(OBJDIR)/lib_discrete_extra.o
	$(CPP) $(CPPFLAGS) $(OFLAGS) -o $@ $^ $(LIBS) -ltiff -ljpeg -lIlmImf -lHalf -lfftw3 -lfftw3f #-lfftw3q

$(OBJDIR)/lib_check_extrema.o : $(SRCDIR)/lib_check_extrema.c $(OBJDIR)/lib_scalespace.o $(OBJDIR)/lib_util.o
	$(CC) $(CFLAGS) $(OFLAGS) -c -o $@ $^


# checking the presence of extrema

$(BINDIR)/grab_ss_neighborhood : $(SRCDIR)/grab_ss_neighborhood.cpp $(OBJ) $(OBJDIR)/iio.o $(OBJDIR)/io_png.o  $(OBJDIR)/io_exr.o $(OBJDIR)/lib_check_extrema.o $(OBJDIR)/lib_fourier.o $(OBJDIR)/lib_dense_anatomy.o $(OBJDIR)/lib_discrete_extra.o
	$(CPP) $(CPPFLAGS) $(OFLAGS) -o $@ $^ $(LIBS) -ltiff -ljpeg -lIlmImf -lHalf -lfftw3 -lfftw3f #-lfftw3q


$(BINDIR)/check_extremum_in_cube : $(SRCDIR)/check_extremum_in_cube.cpp $(OBJ) $(OBJDIR)/iio.o $(OBJDIR)/io_png.o  $(OBJDIR)/io_exr.o  $(OBJDIR)/lib_check_extrema.o  $(OBJDIR)/lib_fourier.o $(OBJDIR)/lib_dense_anatomy.o $(OBJDIR)/lib_discrete_extra.o
	$(CPP) $(CPPFLAGS) $(OFLAGS) -o $@ $^ $(LIBS) -ltiff -ljpeg -lIlmImf -lHalf -lfftw3 -lfftw3f #-lfftw3q

$(BINDIR)/check_extrema_in_scalespace : $(SRCDIR)/check_extrema_in_scalespace.cpp $(OBJ) $(OBJDIR)/iio.o $(OBJDIR)/io_png.o  $(OBJDIR)/io_exr.o  $(OBJDIR)/lib_check_extrema.o  $(OBJDIR)/lib_fourier.o $(OBJDIR)/lib_dense_anatomy.o $(OBJDIR)/lib_discrete_extra.o
	$(CPP) $(CPPFLAGS) $(OFLAGS) -o $@ $^ $(LIBS) -ltiff -ljpeg -lIlmImf -lHalf -lfftw3 -lfftw3f #-lfftw3q

$(BINDIR)/check_extrema_from_scalespace : $(SRCDIR)/check_extrema_from_scalespace.cpp $(OBJ) $(OBJDIR)/iio.o $(OBJDIR)/io_png.o  $(OBJDIR)/io_exr.o  $(OBJDIR)/lib_check_extrema.o  $(OBJDIR)/lib_fourier.o $(OBJDIR)/lib_dense_anatomy.o $(OBJDIR)/lib_discrete_extra.o
	$(CPP) $(CPPFLAGS) $(OFLAGS) -o $@ $^ $(LIBS) -ltiff -ljpeg -lIlmImf -lHalf -lfftw3 -lfftw3f #-lfftw3q

$(BINDIR)/descriptors_from_scalespace : $(SRCDIR)/descriptors_from_scalespace.cpp $(OBJ) $(OBJDIR)/iio.o $(OBJDIR)/io_png.o  $(OBJDIR)/io_exr.o  $(OBJDIR)/lib_check_extrema.o  $(OBJDIR)/lib_fourier.o $(OBJDIR)/lib_dense_anatomy.o $(OBJDIR)/lib_discrete_extra.o
	$(CPP) $(CPPFLAGS) $(OFLAGS) -o $@ $^ $(LIBS) -ltiff -ljpeg -lIlmImf -lHalf -lfftw3 -lfftw3f #-lfftw3q

$(BINDIR)/detection_from_scalespace : $(SRCDIR)/detection_from_scalespace.c $(OBJ) $(OBJDIR)/iio.o $(OBJDIR)/io_png.o  $(OBJDIR)/lib_fourier.o $(OBJDIR)/lib_dense_anatomy.o $(OBJDIR)/lib_discrete_extra.o
	$(CC) $(CFLAGS) $(OFLAGS) -o $@ $^ $(LIBS) -ltiff -ljpeg -lfftw3 -lfftw3f #-lfftw3q


#-------------------------------------------------------------------------------------
#  clean
#
cleanobj:
	-rm -f $(OBJ)

clean: cleanobj
	-rm -f $(BIN)
