TEMPLATE = app
CONFIG += console c++14
CONFIG -= app_bundle
CONFIG -= qt

INCLUDEPATH += $$PWD \
    $$PWD/src/ \
    $$PWD/src/Flux/ \
    $$PWD/src/Problem/ \
    $$PWD/src/Integrator/ \
    $$PWD/src/Mesh/ \
    $$PWD/src/Solver/

SOURCES += RKDG2D.cpp \
    src/Mesh/Cell.cpp \
    src/Mesh/Edge.cpp \
    src/Mesh/EdgeInternal.cpp \
    src/Mesh/Mesh2D.cpp \
    src/Mesh/EdgeBoundary.cpp \
    src/Mesh/EdgeBoundaryInfty.cpp \
    src/Problem/Problem.cpp \
    src/Solver/Solver.cpp \
    src/Flux/Flux.cpp \
    src/Flux/FluxLLF.cpp

HEADERS += src/numvector.h \
    src/Mesh/Mesh2D.h \
    src/Mesh/Point.h \
    src/Mesh/Edge.h \
    src/Mesh/EdgeInternal.h \
    src/Mesh/Cell.h \
    src/Mesh/EdgeBoundary.h \
    src/Mesh/EdgeBoundaryInfty.h \
    src/Problem/Problem.h \
    src/Solver/Solver.h \
    src/Flux/Flux.h \
    src/Flux/FluxLLF.h
