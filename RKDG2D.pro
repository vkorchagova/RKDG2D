TEMPLATE = app
CONFIG += console c++11
CONFIG -= app_bundle
CONFIG -= qt

INCLUDEPATH += $$PWD \
    $$PWD/src/ \
    $$PWD/src/Flux/ \
    $$PWD/src/Problem/ \
    $$PWD/src/Integrator/ \
    $$PWD/src/Mesh/

SOURCES += RKDG2D.cpp \
    #src/Flux/Flux.cpp \
    #src/Flux/FluxLLF.cpp \
    src/Mesh/Edge.cpp \
    src/Mesh/Mesh2D.cpp \
    src/Mesh/Cell.cpp \
    src/Problem/Problem.cpp

HEADERS += src/numvector.h \
    #src/Flux/Flux.h \
    #src/Flux/FluxLLF.h \
    #src/Integrator/gaussintegrator.h \
    src/Mesh/Edge.h \
    src/Mesh/Mesh2D.h \
    src/numvector.h \
    src/Mesh/Cell.h \
    src/Mesh/Point.h \
    src/Problem/Problem.h
