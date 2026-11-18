/* version.h: library version (CMakeLists.txt reads C2FXP_VERSION from here).
 * CODEC2_VERSION_* give the upstream libcodec2 API level this library implements (codec2.h of Codec2 1.2),
 * so programs that check them build unchanged. C2FXP_VERSION is this library's own version. */
#ifndef CODEC2_HAVE_VERSION
#define CODEC2_HAVE_VERSION

#define CODEC2_VERSION_MAJOR 1
#define CODEC2_VERSION_MINOR 2
#define CODEC2_VERSION_PATCH 0
#define CODEC2_VERSION "1.2.0"

#define CODEC2_FIXED_POINT 1
#define C2FXP_VERSION "0.1.0"

#endif
