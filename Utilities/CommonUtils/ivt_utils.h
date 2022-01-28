#ifndef __IVT_UTILS_H__
#define __IVT_UTILS_H__

#ifdef __IVT_FOUND__

#include <ittnotify.h>
#define IVT_PAUSE __itt_pause()
#define IVT_RESUME __itt_resume()

#else

#define IVT_PAUSE
#define IVT_RESUME

#endif

#endif