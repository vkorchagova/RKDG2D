TEMPLATE = app
CONFIG += console c++11
CONFIG -= app_bundle
CONFIG -= qt

SOURCES += RKDG2D.cpp \
    Mesh2D.cpp \
    Flux.cpp \
    Problem.cpp \
    gaussintegrator.cpp

HEADERS += \
    Mesh2D.h \
    numvector.h \
    Flux.h \
    Problem.h \
    gaussintegrator.h
