/*
MIT License

Copyright (c) 2025 hls6432

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
*/

#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdlib>
#include <tuple>
#include <vector>
#include <time.h>
#include <omp.h>

using namespace std;

// Universal constants/fixed values
const double k_B = 0.0138;
const double k_e = 8990;
const double mass[5] = {66.3, 38.2, 58.8, 26.6, 1.67};
const double charge[5] = {0, 0.160, -0.160, -0.1312, 0.0656};
const double sigma[5] = {0.340, 0.216, 0.483, 0.317, 0};
const double epsilon[5] = {1.67, 2.45, 0.0889, 1.08, 0};
const string solvnames[2] = {"Argon", "Water"};
const double solvmasses[2] = {mass[0], mass[3] + 2 * mass[4]};
const double latconst = 0.564;
const double rOH = 0.100;
const double thetaHOH = 1.91;
const double rHH = 2 * rOH * sin(0.5 * thetaHOH);
const double watera = 933200;
const double waterb = 228500;
const double waterc = -147100;
const double waterd = 77700;

// Simulation parameters
const double dt = 0.0005;
const double rc = 0.9;
const double rc2 = rc * rc;
const int iterations = 500000;
const int interval = 500;

// Setup parameters
const int solvtype = 1;
const double tempinit = 323;
const int nsaltlen = 3;
const double bsizes[2] = {29.6, 2.60};
const double soluspace = 0.5 * latconst;
const double solvspaces[2] = {3.34, 0.300};

const string solvname = solvnames[solvtype];
const double bsize = bsizes[solvtype];
const double solvspace = solvspaces[solvtype];
const double solvmass = solvmasses[solvtype];
const int nsalt = nsaltlen * nsaltlen * nsaltlen;
int nsolvlen = int(bsize / solvspace);
int nsolv, npart;

// Simulation variables
vector<array<double, 3>> pos = {};
vector<array<double, 3>> ppos = {};
vector<array<double, 3>> vel = {};
vector<array<double, 3>> force = {};
vector<int> type = {};
vector<int> waterO = {}, waterH1 = {}, waterH2 = {};

double temp = tempinit;
double en = 0;

// Output
string datafile = "data.txt";
ofstream pdata;

// Uniformly samples a random number from 0 to 1
double random() {
    return ((double) rand()) / RAND_MAX;
}

// Generates a random orientation for a water molecule and returns the displacement vectors for the two O-H bonds
tuple<double, double, double, double, double, double> randwater() {
    double u1 = random(), u2 = random(), u3 = random();
    double a = sqrt(1 - u1), b = sqrt(u1);
    double w = a * sin(2 * M_PI * u2), x = a * cos(2 * M_PI * u2), y = b * sin(2 * M_PI * u3), z = b * cos(2 * M_PI * u3);
    double a11 = 1 - 2 * (y * y + z * z), a12 = 2 * (x * y - z * w); // a13 = 2 * (x * z + y * w);
    double a21 = 2 * (x * y + z * w), a22 = 1 - 2 * (x * x + z * z); // a23 = 2 * (y * z - x * w);
    double a31 = 2 * (x * z - y * w), a32 = 2 * (y * z + x * w); // a33 = 1 - 2 * (x * x + y * y);
    double x1 = a11 * rOH, y1 = a21 * rOH, z1 = a31 * rOH;
    double x2 = a11 * rOH * cos(thetaHOH) + a12 * rOH * sin(thetaHOH);
    double y2 = a21 * rOH * cos(thetaHOH) + a22 * rOH * sin(thetaHOH);
    double z2 = a31 * rOH * cos(thetaHOH) + a32 * rOH * sin(thetaHOH);
    return make_tuple(x1, y1, z1, x2, y2, z2);
}

// Initialisation
void init() {
    pos.clear(); ppos.clear(); vel.clear(); force.clear();
    type.clear(); waterO.clear(); waterH1.clear(); waterH2.clear();

    // Initialise solvent molecules
    double meanvx = 0, meanvy = 0, meanvz = 0, meanv2 = 0;
    double lim = soluspace * (nsaltlen + 1) * 0.5;
    nsolv = 0;
    for (int i = 0; i < nsolvlen; i++) {
        for (int j = 0; j < nsolvlen; j++) {
            for (int k = 0; k < nsolvlen; k++) {
                double x = solvspace * (i - (nsolvlen - 1) * 0.5);
                double y = solvspace * (j - (nsolvlen - 1) * 0.5);
                double z = solvspace * (k - (nsolvlen - 1) * 0.5);
                if (abs(x) > lim || abs(y) > lim || abs(z) > lim) {
                    x += 0.5 * bsize; y += 0.5 * bsize; z += 0.5 * bsize;
                    double vx = random() - 0.5, vy = random() - 0.5, vz = random() - 0.5;
                    double px = x - vx * dt, py = y - vy * dt, pz = z - vz * dt;
                    meanvx += vx; meanvy += vy; meanvz += vz;
                    meanv2 += vx * vx + vy * vy + vz * vz;
                    if (solvtype == 0) {
                        array<double, 3> p = {x, y, z}, v = {vx, vy, vz};
                        pos.push_back(p); vel.push_back(v); type.push_back(0);
                    } else if (solvtype == 1) {
                        double dx1, dy1, dz1, dx2, dy2, dz2;
                        tie(dx1, dy1, dz1, dx2, dy2, dz2) = randwater();
                        double x1 = x + dx1, y1 = y + dy1, z1 = z + dz1;
                        double x2 = x + dx2, y2 = y + dy2, z2 = z + dz2;
                        array<double, 3> p = {x, y, z}, p1 = {x1, y1, z1}, p2 = {x2, y2, z2};
                        array<double, 3> v = {vx, vy, vz};
                        pos.push_back(p); pos.push_back(p1); pos.push_back(p2);
                        vel.push_back(v); vel.push_back(v); vel.push_back(v);
                        type.push_back(3); type.push_back(4); type.push_back(4);
                        waterO.push_back(3*nsolv); waterH1.push_back(3*nsolv+1); waterH2.push_back(3*nsolv+2);
                    }
                    nsolv++;
                }
            }
        }
    }
    meanvx /= nsolv; meanvy /= nsolv; meanvz /= nsolv; meanv2 /= nsolv;
    double scale = sqrt(3 * tempinit * k_B / (solvmass * meanv2));
    if (solvtype == 0) {
        for (int i = 0; i < nsolv; i++) {
            array<double, 3> v = vel.at(i);
            double newvx = (v[0] - meanvx) * scale, newvy = (v[1] - meanvy) * scale, newvz = (v[2] - meanvz) * scale;
            array<double, 3> newv = {newvx, newvy, newvz};
            vel.at(i) = newv;
            array<double, 3> p = pos.at(i);
            array<double, 3> pp = {p[0] - newvx * dt, p[1] - newvy * dt, p[2] - newvz * dt};
            ppos.push_back(pp);
        }
    } else if (solvtype == 1) {
        for (int i = 0; i < nsolv; i++) {
            array<double, 3> v = vel.at(3*i);
            double newvx = (v[0] - meanvx) * scale, newvy = (v[1] - meanvy) * scale, newvz = (v[2] - meanvz) * scale;
            array<double, 3> newv = {newvx, newvy, newvz};
            for (int j = 3*i; j < 3*i+3; j++) {
                vel.at(j) = newv;
                array<double, 3> p = pos.at(j);
                array<double, 3> pp = {p[0] - newvx * dt, p[1] - newvy * dt, p[2] - newvz * dt};
                ppos.push_back(pp);
            }
        }
    }

    // Initialise salt molecules
    for (int i = 0; i < nsaltlen; i++) {
        for (int j = 0; j < nsaltlen; j++) {
            for (int k = 0; k < nsaltlen; k++) {
                double x = soluspace * (i - (nsaltlen - 1) * 0.5) + 0.5 * bsize;
                double y = soluspace * (j - (nsaltlen - 1) * 0.5) + 0.5 * bsize;
                double z = soluspace * (k - (nsaltlen - 1) * 0.5) + 0.5 * bsize;
                array<double, 3> p = {x, y, z}, v = {0, 0, 0}, pp = {x, y, z};
                pos.push_back(p);
                vel.push_back(v);
                ppos.push_back(pp);
                type.push_back((i + j + k) % 2 + 1);
            }
        }
    }
    npart = pos.size();
    array<double, 3> f = {0, 0, 0};
    for (int i = 0; i < npart; i++)
        force.push_back(f);
}

// Calculate forces, parallel processed using OpenMP
void calcforces() {
    en = 0;
    #pragma omp parallel for
    for (int i = 0; i < npart; i++)
        force.at(i) = {0, 0, 0};
    #pragma omp parallel
    {
        vector<array<double, 3>> force_t(npart, {0, 0, 0});
        double en_t = 0;

        // Calculate Lennard-Jones forces and Coloumb forces
        #pragma omp for schedule(dynamic)
        for (int i = 0; i < npart - 1; i++) {
            for (int j = i + 1; j < npart; j++) {
                int typei = type.at(i), typej = type.at(j);
                array<double, 3> pi = pos.at(i), pj = pos.at(j);
                double sig = (sigma[typei] + sigma[typej]) * 0.5;
                double eps = sqrt(epsilon[typei] * epsilon[typej]);
                double rx = pi[0] - pj[0]; double ry = pi[1] - pj[1]; double rz = pi[2] - pj[2];
                rx = rx - bsize * round(rx / bsize); ry = ry - bsize * round(ry / bsize); rz = rz - bsize * round(rz / bsize);
                double r2 = rx * rx + ry * ry + rz * rz;
                if (0 < r2 && r2 < rc2) {
                    double r2i = (sig * sig) / r2;
                    double r6i = r2i * r2i * r2i;
                    double f = 48 * eps * r6i * (r6i - 0.5) / r2;
                    force_t.at(i)[0] += f * rx; force_t.at(i)[1] += f * ry; force_t.at(i)[2] += f * rz;
                    force_t.at(j)[0] -= f * rx; force_t.at(j)[1] -= f * ry; force_t.at(j)[2] -= f * rz;
                    double rc2i = (sig * sig) / rc2;
                    double rc6i = rc2i * rc2i * rc2i;
                    en_t += 4 * eps * (r6i * (r6i - 1) - rc6i * (rc6i - 1));
                }
                if (charge[typei] != 0 && charge[typej] != 0) {
                    double r = sqrt(r2);
                    double f = k_e * charge[typei] * charge[typej] / (r * r2);
                    force_t.at(i)[0] += f * rx; force_t.at(i)[1] += f * ry; force_t.at(i)[2] += f * rz;
                    force_t.at(j)[0] -= f * rx; force_t.at(j)[1] -= f * ry; force_t.at(j)[2] -= f * rz;
                    en_t += k_e * charge[typei] * charge[typej] / r;
                }
            }
        }
        #pragma omp critical
        {
            for (int i = 0; i < npart; i++) {
                force.at(i)[0] += force_t.at(i)[0];
                force.at(i)[1] += force_t.at(i)[1];
                force.at(i)[2] += force_t.at(i)[2];
            }
            en += en_t;
        }
    }

    // Calculate intramolecular forces between water molecule atoms
    if (solvtype == 1) {
        for (int i = 0; i < nsolv; i++) {
            int oi = waterO.at(i); int hi = waterH1.at(i); int hj = waterH2.at(i);
            array<double, 3> poi = pos[oi], phi = pos[hi], phj = pos[hj];
            double x1 = poi[0] - phi[0], y1 = poi[1] - phi[1], z1 = poi[2] - phi[2];
            double x2 = poi[0] - phj[0], y2 = poi[1] - phj[1], z2 = poi[2] - phj[2];
            double x3 = phi[0] - phj[0], y3 = phi[1] - phj[1], z3 = phi[2] - phj[2];
            double r1 = sqrt(x1 * x1 + y1 * y1 + z1 * z1) - rOH;
            double r2 = sqrt(x2 * x2 + y2 * y2 + z2 * z2) - rOH;
            double r3 = sqrt(x3 * x3 + y3 * y3 + z3 * z3) - rHH;
            double f1 = (watera * r1 + waterc * r3 + waterd * r2) / (r1 + rOH);
            double f2 = (watera * r2 + waterc * r3 + waterd * r1) / (r2 + rOH);
            double f3 = (waterb * r3 + waterc * (r1 + r2)) / (r3 + rHH);
            force.at(oi)[0] -= f1 * x1 + f2 * x2; force.at(oi)[1] -= f1 * y1 + f2 * y2; force.at(oi)[2] -= f1 * z1 + f2 * z2;
            force.at(hi)[0] -= -f1 * x1 + f3 * x3; force.at(hi)[1] -= -f1 * y1 + f3 * y3; force.at(hi)[2] -= -f1 * z1 + f3 * z3;
            force.at(hj)[0] -= -f2 * x2 - f3 * x3; force.at(hj)[1] -= -f2 * y2 - f3 * y3; force.at(hj)[2] -= -f2 * z2 - f3 * z3;
            en += 0.5 * watera * (r1 * r1 + r2 * r2) + 0.5 * waterb * r3 * r3 + waterc * (r1 + r2) * r3 + waterd * r1 * r2;
        }
    }
}

// Integrate equations of motion
void intmotion() {
    for (int i = 0; i < npart; i++) {
        double m = mass[type.at(i)];
        array<double, 3> p = pos.at(i), pp = ppos.at(i), f = force.at(i);
        double x = p[0], y = p[1], z = p[2];
        double xp = pp[0], yp = pp[1], zp = pp[2];
        // Verlet integration
        double xn = 2 * x - xp + dt * dt * f[0] / m;
        double yn = 2 * y - yp + dt * dt * f[1] / m;
        double zn = 2 * z - zp + dt * dt * f[2] / m;
        double vx = (xn - xp) / (2 * dt);
        double vy = (yn - yp) / (2 * dt);
        double vz = (zn - zp) / (2 * dt);
        en += 0.5 * mass[type[i]] * (vx * vx + vy * vy + vz * vz);
        ppos.at(i)[0] = x; ppos.at(i)[1] = y; ppos.at(i)[2] = z;
        pos.at(i)[0] = xn; pos.at(i)[1] = yn; pos.at(i)[2] = zn;
        vel.at(i)[0] = vx; vel.at(i)[1] = vy; vel.at(i)[2] = vz;
    }

    // Update quantities
    double meanv2 = 0;
    for (int i = 0; i < nsolv; i++) {
        array<double, 3> v = vel.at(i);
        double vx = v[0], vy = v[1], vz = v[2];
        meanv2 += v[0] * v[0] + v[1] * v[1] + v[2] * v[2];
    }
    meanv2 /= nsolv;
    temp = solvmass * meanv2 / (3 * k_B);
}

// Initialise data file
void initdata() {
    pdata.open(datafile);
    pdata << solvname << " " << tempinit << " " << nsaltlen << " " << soluspace << " " << solvspace << "\n";
    pdata << bsize << " " << rc << " " << dt << " " << iterations << " " << interval << "\n";
}

// Log current data into data file
void logdata() {
    for (int i = 0; i < npart; i++) {
        pdata << pos[i][0] << " " << pos[i][1] << " " << pos[i][2] << " ";
        pdata << vel[i][0] << " " << vel[i][1] << " " << vel[i][2] << " ";
        pdata << type[i];
        if (i < npart - 1)
            pdata << ", ";
    }
    pdata << " | " << en << " | " << temp;
    pdata << "\n";
}

// Main loop
int main() {
    init();
    initdata();
    cout << "Max threads: " << omp_get_max_threads() << "\n";
    cout << "Initial temperature: " << tempinit << " K\n";
    cout << "Solvent type: " << solvname << "\n";
    cout << "Number of salt ions: " << nsalt << "\n";
    cout << "Number of solvent particles: " << nsolv << "\n";
    cout << "Total atoms: " << npart << "\n";
    cout << "Box size: " << bsize << " nm\n";
    cout << "Timestep: " << dt << " ps\n";
    cout << "Cutoff radius: " << rc << " nm\n";
    cout << "Iterations: " << iterations << "\n";
    cout << "Sampling interval: " << interval << "\n";
    cout << "\n";
    time_t start = time(0);
    for (int j = 0; j < iterations; j++) {
        calcforces();
        intmotion();
        bool stop = false;
        for (int i = 0; i < npart; i++) {
            if (abs(vel[i][0]) > 1000 || abs(vel[i][1]) > 1000 || abs(vel[i][2]) > 1000) {
                cout << "Velocity has become big, stopping simulation...\n";
                stop = true;
            }
        }
        if (stop)
            break;
        if (j % interval == 0) {
            cout << "Iteration " << j << "/" << iterations << ". ";
            double secs = (double) difftime(time(0), start);
            cout << "Elapsed time: " << round(secs * 100) / 100 << " seconds\n";
            logdata();
        }
    }
    pdata << "end";
    cout << "Complete! Total elapsed time: " << round(((double) difftime(time(0), start)) / 60 * 100) / 100 << "minutes \n";
    pdata.close();
    return 0;
}
