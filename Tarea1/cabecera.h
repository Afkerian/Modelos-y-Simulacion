#ifndef CABECERA_H_INCLUDED
#define CABECERA_H_INCLUDED

#include <iostream>
using namespace std;

class vetor {
public:
    double x, y, z;

    // Constructores
    vetor() : x(0.0), y(0.0), z(0.0) {}
    vetor(double xi, double yi, double zi) : x(xi), y(yi), z(zi) {}

    // Producto punto
    double operator*(const vetor& v) const {
        return x * v.x + y * v.y + z * v.z;
    }

    // Producto por un escalar
    vetor operator*(double scalar) const {
        return vetor(x * scalar, y * scalar, z * scalar);
    }

    // Resta de vectores
    vetor operator-(const vetor& v) const {
        return vetor(x - v.x, y - v.y, z - v.z);
    }

    // Módulo del vector
    double modulo() const {
        return sqrt(x * x + y * y + z * z);
    }

    // Normalizar el vector
    void normalizar() {
        double m = modulo();
        if (m > 0) {
            x /= m;
            y /= m;
            z /= m;
        }
    }

    // Producto cruz
    vetor cruz(const vetor& v) const {
        return vetor(
            y * v.z - z * v.y,
            z * v.x - x * v.z,
            x * v.y - y * v.x
        );
    }
};


class point {
public:
    double x, y, z;

    // Constructores
    point() : x(0.0), y(0.0), z(0.0) {}
    point(double xi, double yi, double zi) : x(xi), y(yi), z(zi) {}

    // Método para inicializar
    void iniciar(double xi = 0.0, double yi = 0.0, double zi = 0.0) {
        x = xi;
        y = yi;
        z = zi;
    }

    // Sobrecarga del operador '-' para obtener un vector entre dos puntos
    vetor operator-(const point& p) const {
        return vetor(x - p.x, y - p.y, z - p.z);
    }

    void printPoint() const {
        cout << "(" << x << ", " << y << ", " << z << ")" << endl;
    }
};



class plane {
public:
    point *p;
    int NP;
    vetor n;

    plane() : p(NULL), NP(0), n() {}

    void iniciar() {
        p = NULL;
        NP = 0;
        n = vetor();
    }

    void NewPoints(int N) {
        point *tp = new point[NP + N];
        for (int P = 0; P < NP; P++) {
            tp[P] = p[P];
        }
        for (int P = NP; P < NP + N; P++) {
            tp[P].iniciar();
        }
        delete[] p;
        p = tp;
        NP += N;
    }

    void calcularNormal() {
        if (NP >= 3) {
            vetor v1 = p[1] - p[0];
            vetor v2 = p[2] - p[0];
            n = v1.cruz(v2);
            n.normalizar();
        } else {
            cout << "No hay suficientes puntos para calcular la normal." << endl;
        }
    }
};


class room {
public:
    plane *p;
    int NP;

    room() : p(NULL), NP(0) {}

    void NewPlane(int N) {
        plane *tp = new plane[NP + N];
        for (int P = 0; P < NP; P++) {
            tp[P] = p[P];
        }
        for (int P = NP; P < NP + N; P++) {
            tp[P].iniciar();
        }
        delete[] p;
        p = tp;
        NP += N;
    }
};

class Simulation {
public:
    Room room;        // El cubo
    Ray ray;          // El rayo
    int reflections;  // Número de reflexiones

    // Constructor
    Simulation(int numReflections) : reflections(numReflections) {
        // Generar dirección aleatoria para el rayo
        srand(time(0));
        Vector v(
            ((double)rand() / RAND_MAX) - 0.5,
            ((double)rand() / RAND_MAX) - 0.5,
            ((double)rand() / RAND_MAX) - 0.5
        );
        v.normalize();

        ray = Ray(Point(0.0, 0.0, 0.0), v);
    }

    // Método para ejecutar la simulación
    void run() {
        for (int i = 0; i < reflections; ++i) {
            double minT = DBL_MAX;
            Point intersectionPoint;
            int planeIndex = -1;

            // Buscar el plano más cercano con el que el rayo intersecta
            for (int j = 0; j < 6; ++j) {
                Point tempIntersection;
                double t;
                if (ray.intersectPlane(room.planes[j], tempIntersection, t)) {
                    if (t < minT) {
                        minT = t;
                        intersectionPoint = tempIntersection;
                        planeIndex = j;
                    }
                }
            }

            if (planeIndex == -1) {
                std::cout << "El rayo ha salido del cubo." << std::endl;
                break;
            }

            // Actualizar la posición y dirección del rayo
            ray.origin = intersectionPoint;
            ray.reflect(room.planes[planeIndex].normal);

            // Imprimir la posición actual
            std::cout << "Reflexión " << i + 1 << ": ";
            ray.origin.print();
            std::cout << std::endl;
        }
    }
};

class Ray {
public:
    Point origin;     // Origen del rayo
    Vector direction; // Dirección del rayo

    // Constructores
    Ray() {}
    Ray(const Point& o, const Vector& d) : origin(o), direction(d) {}

    // Métodos
    bool intersectPlane(const Plane& plane, Point& intersection, double& t) const {
        double denom = plane.normal.dot(direction);
        if (fabs(denom) > 1e-6) {
            Vector p0l0 = plane.p1 - origin;
            t = plane.normal.dot(p0l0) / denom;
            if (t >= 1e-6) {
                intersection = Point(
                    origin.x + direction.x * t,
                    origin.y + direction.y * t,
                    origin.z + direction.z * t
                );
                return true;
            }
        }
        return false;
    }

    void reflect(const Vector& normal) {
        direction = direction - normal * (2 * direction.dot(normal));
        direction.normalize();
    }
};


#endif // CABECERA_H_INCLUDED



#endif // CABECERA_H_INCLUDED
