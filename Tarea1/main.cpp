// main.cpp

#include "cabecera.h"
#include <iostream>
#include <cstdlib>
#include <ctime>
#include <float.h>

int main() {
    // Crear una simulación con 50 reflexiones
    Simulation sim(50);

    // Ejecutar la simulación
    sim.run();

    return 0;
}
