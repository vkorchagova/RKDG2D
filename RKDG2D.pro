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
    $$PWD/src/Indicator/ \
    $$PWD/src/Limiter/

SOURCES += RKDG2D.cpp \
    src/Mesh/Cell.cpp \
    src/Mesh/Edge.cpp \
    src/Mesh/EdgeInternal.cpp \
    src/Mesh/Mesh2D.cpp \
    src/Mesh/EdgeBoundary.cpp \
    src/Mesh/EdgeBoundaryInfty.cpp \
    src/Mesh/EdgeBoundaryOpen.cpp \
    src/Mesh/EdgeBoundarySlip.cpp \
    src/Problem/Problem.cpp \
    src/Solver/Solver.cpp \
    src/Flux/Flux.cpp \
    src/Flux/FluxLLF.cpp \
    src/Flux/FluxHLL.cpp \
    src/Flux/FluxHLLC.cpp \
    src/defs/defs.cpp \
    src/Indicator/Indicator.cpp \
    src/Indicator/IndicatorEverywhere.cpp \
    src/Indicator/IndicatorNowhere.cpp \
    src/Indicator/IndicatorKXRCF.cpp \
    src/Limiter/Limiter.cpp \
    src/Limiter/LimiterFinDiff.cpp \
    src/Limiter/LimiterMUSCL.cpp \
    src/Limiter/LimiterWENOS.cpp


HEADERS += src/numvector.h \
    src/Mesh/Mesh2D.h \
    src/Mesh/Point.h \
    src/Mesh/Edge.h \
    src/Mesh/EdgeInternal.h \
    src/Mesh/Cell.h \
    src/Mesh/EdgeBoundary.h \
    src/Mesh/EdgeBoundaryInfty.h \
    src/Mesh/EdgeBoundaryOpen.h \
    src/Mesh/EdgeBoundarySlip.h \
    src/Problem/Problem.h \
    src/Solver/Solver.h \
    src/Flux/Flux.h \
    src/Flux/FluxLLF.h \
    src/Flux/FluxHLL.h \
    src/Flux/FluxHLLC.h \
    src/defs/defs.h \
    src/Indicator/Indicator.h \
    src/Indicator/IndicatorEverywhere.h \
    src/Indicator/IndicatorNowhere.h \
    src/Indicator/IndicatorKXRCF.h \
    src/Limiter/Limiter.h \
    src/Limiter/LimiterFinDiff.h \
    src/Limiter/LimiterMUSCL.h \
    src/Limiter/LimiterWENOS.h
