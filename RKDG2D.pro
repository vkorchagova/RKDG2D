TEMPLATE = app
CONFIG += console c++14
CONFIG -= app_bundle
CONFIG -= qt

INCLUDEPATH += $$PWD \
    $$PWD/src/ \
    $$PWD/src/defs/ \
    $$PWD/src/Flux/ \
    $$PWD/src/Problem/ \
    $$PWD/src/Integrator/ \
    $$PWD/src/Mesh/ \
    $$PWD/src/Solver/ \
    $$PWD/src/Indicator/

SOURCES += RKDG2D.cpp \
    src/Mesh/Cell.cpp \
    src/Mesh/Edge.cpp \
    src/Mesh/EdgeInternal.cpp \
    src/Mesh/Mesh2D.cpp \
    src/Mesh/EdgeBoundary.cpp \
    src/Mesh/EdgeBoundaryInfty.cpp \
    src/Mesh/EdgeBoundarySlip.cpp \
    src/Problem/Problem.cpp \
    src/Solver/Solver.cpp \
    src/Flux/Flux.cpp \
    src/Flux/FluxLLF.cpp \
    src/Flux/FluxHLL.cpp \
    src/defs/defs.cpp \
    src/Indicator/Indicator.cpp \
    src/Indicator/IndicatorKXRCF.cpp


HEADERS += src/numvector.h \
    src/Mesh/Mesh2D.h \
    src/Mesh/Point.h \
    src/Mesh/Edge.h \
    src/Mesh/EdgeInternal.h \
    src/Mesh/Cell.h \
    src/Mesh/EdgeBoundary.h \
    src/Mesh/EdgeBoundaryInfty.h \
    src/Mesh/EdgeBoundarySlip.h \
    src/Problem/Problem.h \
    src/Solver/Solver.h \
    src/Flux/Flux.h \
    src/Flux/FluxLLF.h \
    src/Flux/FluxHLL.h \
    src/defs/defs.h \
    src/Indicator/Indicator.h \
    src/Indicator/IndicatorKXRCF.h
