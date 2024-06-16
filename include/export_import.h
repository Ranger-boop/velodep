#ifndef EXPORT_IMPORT_H_
#define EXPORT_IMPORT_H_

#ifdef _WIN32
    #ifdef BUILDING_DLL
        #define API __declspec(dllexport)
    #else
        #define API __declspec(dllimport)
    #endif
#else
    #define API
#endif

#endif  // EXPORT_IMPORT_H_