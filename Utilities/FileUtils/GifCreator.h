// Modified from gif.h by Charlie Tangora 
// Github of original file: https://github.com/charlietangora/gif-h
//
// This file offers a simple, very limited way to create animated GIFs directly in code.
//
// Those looking for particular cleverness are likely to be disappointed; it's pretty
// much a straight-ahead implementation of the GIF format with optional Floyd-Steinberg
// dithering. (It does at least use delta encoding - only the changed portions of each
// frame are saved.)
//
// So resulting files are often quite large. The hope is that it will be handy nonetheless
// as a quick and easily-integrated way for programs to spit out animations.
//
// Only RGBA8 is currently supported as an input format. (The alpha is ignored.)
//
// If capturing a buffer with a bottom-left origin (such as OpenGL), define GIF_FLIP_VERT
// to automatically flip the buffer data when writing the image (the buffer itself is
// unchanged.
//
// USAGE:
// 1. Create a GifWriter struct.
// 2. Pass it to GifBegin() to initialize and write the header.
// 3. Pass subsequent frames to GifWriteFrame().
// 4. Call GifEnd() to close the file handle and free memory.
//

#ifndef _Gif_Creator_h_
#define _Gif_Creator_h_

namespace GifCreator
{
    struct GifWriter
    {
        FILE* f;
        uint8_t* oldImage;
        bool firstFrame;
    };

    // Creates a gif file.
    // The input GIFWriter is assumed to be uninitialized.
    bool GifBegin(
        GifWriter* writer,
        const char* filename,
        uint32_t width, uint32_t height,
        bool is_repeatable, // whether the gif is played repeatedly
        int32_t bitDepth = 8,
        bool dither = false
        );

    // Writes out a new frame to a GIF in progress.
    // The delay value is the time between frames in hundredths of a second
    // The GIFWriter should have been created by GIFBegin.
    // AFAIK, it is legal to use different bit depths for different frames of an image -
    // this may be handy to save bits in animations that don't change much.
    bool GifWriteFrame(
        GifWriter* writer,
        const uint8_t* image,
        uint32_t width,
        uint32_t height,
        uint32_t delay,
        int bitDepth = 8,
        bool dither = false
        );

    // Writes the EOF code, closes the file handle, and frees temp memory used by a GIF.
    // Many if not most viewers will still display a GIF properly if the EOF code is missing,
    // but it's still a good idea to write it out.
    bool GifEnd(GifWriter* writer);
}

#endif