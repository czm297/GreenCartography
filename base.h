#pragma once


#ifndef MAPSTYLETRANSFORM_EXPORT
#if defined(_MSC_VER) || defined(__CYGWIN__) || defined(__MINGW32__) || defined( __BCPLUSPLUS__)  || defined( __MWERKS__)
# define MAPSTYLETRANSFORM_EXPORT   __declspec(dllexport)
#else
# define MAPSTYLETRANSFORM_EXPORT
#endif
#endif // !MAPSTYLETRANSFORM_EXPORT