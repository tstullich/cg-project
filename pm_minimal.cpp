
/* -------------------------------------------------------------------------
  Minimal (unoptimized) example of PatchMatch. Requires that ImageMagick be installed.

  To improve generality you can:
   - Use whichever distance function you want in dist(), e.g. compare SIFT descriptors computed densely.
   - Search over a larger search space, such as rotating+scaling patches (see MATLAB mex for examples of both)

  To improve speed you can:
   - Turn on optimizations (/Ox /Oi /Oy /fp:fast or -O6 -s -ffast-math -fomit-frame-pointer -fstrength-reduce -msse2 -funroll-loops)
   - Use the MATLAB mex which is already tuned for speed
   - Use multiple cores, tiling the input. See our publication "The Generalized PatchMatch Correspondence Algorithm"
   - Tune the distance computation: manually unroll loops for each patch size, use SSE instructions (see readme)
   - Precompute random search samples (to avoid using rand, and mod)
   - Move to the GPU
  -------------------------------------------------------------------------- */

#include <algorithm>
#include <climits>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <istream>
#include <limits>
#include <random>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unordered_map>
#include <utility>
#include <vector>

#ifndef MAX
#define MAX(a, b) ((a)>(b)?(a):(b))
#define MIN(a, b) ((a)<(b)?(a):(b))
#endif

/* -------------------------------------------------------------------------
   BITMAP: Minimal image class
   ------------------------------------------------------------------------- */

class BITMAP {
public:
    int w, h;
    int *data;

    BITMAP(int w_, int h_) : w(w_), h(h_) { data = new int[w * h]; }

    ~BITMAP() { delete[] data; }

    int *operator[](int y) { return &data[y * w]; }
};

struct pair_hash {
    template<class T1, class T2>
    std::size_t operator()(const std::pair<T1, T2> &p) const {
        auto h1 = std::hash<T1>{}(p.first);
        auto h2 = std::hash<T2>{}(p.second);

        // Mainly for demonstration purposes, i.e. works but is overly simple
        // In the real world, use sth. like boost.hash_combine
        return h1 ^ h2;
    }
};

class Histogram {
public:
    explicit Histogram(const std::vector<int> &samples) {
        hist = std::vector<float>(256);
        for (auto sample : samples) {
            // Decode RGB values of each sample
            int r = sample & 255;
            int g = (sample >> 8) & 255;
            int b = sample >> 16;
            int l = intensity(r, g, b);
            hist[l] += 1.0f;
        }
        // Normalize buckets
        for (float &i : hist) {
            i /= samples.size();
        }
    }

    ~Histogram() {
        std::cout << "Deleting histogram" << std::endl;
        for (auto item : hist) {
            std::cout << item << ", ";
        }
        std::cout << std::endl;
    }

    float operator[](int idx) const { return hist[idx]; }

    int size() const { return hist.size(); }

private:
    // Calculate the color intensity/luminance value
    // Formula taken from https://www.mathworks.com/help/matlab/ref/rgb2gray.html
    int intensity(int r, int g, int b) {
        return static_cast<int>(0.2989f * r + 0.5870f * g + 0.1140f * b);
    }

    std::vector<float> hist;
};

class C2C {
public:
    C2C() {};
    std::unordered_map<std::pair<int, int>, std::vector<int>, pair_hash> data;
};

void check_im() {
    if (system("identify > null.txt") != 0) {
        fprintf(stderr, "ImageMagick must be installed, and 'convert' and 'identify' must be in the path\n");
        exit(1);
    }
}

BITMAP *load_bitmap(const char *filename) {
    //check_im();
    char rawname[256], txtname[256];
    strcpy(rawname, filename);
    strcpy(txtname, filename);
    if (!strstr(rawname, ".")) {
        fprintf(stderr, "Error reading image '%s': no extension found\n", filename);
        exit(1);
    }
    sprintf(strstr(rawname, "."), ".raw");
    sprintf(strstr(txtname, "."), ".txt");
    char buf[256];
    sprintf(buf, "convert %s rgba:%s", filename, rawname);
    if (system(buf) != 0) {
        fprintf(stderr, "Error reading image '%s': ImageMagick convert gave an error\n", filename);
        exit(1);
    }
    sprintf(buf, "identify -format \"%%w %%h\" %s > %s", filename, txtname);
    if (system(buf) != 0) {
        fprintf(stderr, "Error reading image '%s': ImageMagick identify gave an error\n", filename);
        exit(1);
    }
    FILE *f = fopen(txtname, "rt");
    if (!f) {
        fprintf(stderr, "Error reading image '%s': could not read output of ImageMagick identify\n", filename);
        exit(1);
    }
    int w = 0, h = 0;
    if (fscanf(f, "%d %d", &w, &h) != 2) {
        fprintf(stderr, "Error reading image '%s': could not get size from ImageMagick identify\n", filename);
        exit(1);
    }
    fclose(f);
    f = fopen(rawname, "rb");
    BITMAP *ans = new BITMAP(w, h);
    unsigned char *p = (unsigned char *) ans->data;
    for (int i = 0; i < w * h * 4; i++) {
        int ch = fgetc(f);
        if (ch == EOF) {
            fprintf(stderr, "Error reading image '%s': raw file is smaller than expected size %dx%dx4\n", filename, w,
                    h, 4);
            exit(1);
        }
        *p++ = ch;
    }
    fclose(f);
    return ans;
}

int encode_rgb(unsigned char r, unsigned char g, unsigned char b) {
    return ((((b << 16) | 0) | (g << 8)) | r);
}

C2C *load_dump(const char *filename) {
    std::ifstream dump(filename, std::ifstream::binary);
    dump.seekg(0, dump.end);
    int num_floats = dump.tellg() / sizeof(float);
    dump.seekg(0, dump.beg);
    std::vector<float> buf(num_floats);
    dump.read(reinterpret_cast<char *>(buf.data()), buf.size() * sizeof(float));
    auto ans = new C2C();
    for (unsigned int i = 0; i < buf.size(); i += 6) {
        int x = static_cast<int>(floor(buf[i]));
        int y = static_cast<int>(floor(buf[i + 1]));

        float gamma = 1.0f / 2.2f;
        unsigned char r = static_cast<unsigned char>(std::pow(std::clamp(buf[i + 2], 0.0f, 1.0f), gamma) * 255);
        unsigned char g = static_cast<unsigned char>(std::pow(std::clamp(buf[i + 3], 0.0f, 1.0f), gamma) * 255);
        unsigned char b = static_cast<unsigned char>(std::pow(std::clamp(buf[i + 4], 0.0f, 1.0f), gamma) * 255);
        auto xy = std::make_pair(x, y);
        if (ans->data.find(xy) == ans->data.end()) {
            // Insert new pair encoding RGB values in int
            auto pixelbucket = std::vector<int>();
            pixelbucket.emplace_back(encode_rgb(r, g, b));
            ans->data[xy] = pixelbucket;
        } else {
            // Update existing bucket
            ans->data.at(xy).emplace_back(encode_rgb(r, g, b));
        }
    }
    return ans;
}

void save_bitmap(BITMAP *bmp, const char *filename) {
    //check_im();
    char rawname[256];
    strcpy(rawname, filename);
    if (!strstr(rawname, ".")) {
        fprintf(stderr, "Error writing image '%s': no extension found\n", filename);
        exit(1);
    }
    sprintf(strstr(rawname, "."), ".raw");
    char buf[256];
    FILE *f = fopen(rawname, "wb");
    if (!f) {
        fprintf(stderr, "Error writing image '%s': could not open raw temporary file\n", filename);
        exit(1);
    }
    unsigned char *p = (unsigned char *) bmp->data;
    for (int i = 0; i < bmp->w * bmp->h * 4; i++) {
        fputc(*p++, f);
    }
    fclose(f);
    sprintf(buf, "convert -size %dx%d -depth 8 rgba:%s %s", bmp->w, bmp->h, rawname, filename);
    if (system(buf) != 0) {
        fprintf(stderr, "Error writing image '%s': ImageMagick convert gave an error\n", filename);
        exit(1);
    }
}

/* -------------------------------------------------------------------------
   PatchMatch, using L2 distance between upright patches that translate only
   ------------------------------------------------------------------------- */

int patch_w = 7;
int pm_iters = 1;
int rs_max = INT_MAX;

#define XY_TO_INT(x, y) (((y)<<12)|(x))
#define INT_TO_X(v) ((v)&((1<<12)-1))
#define INT_TO_Y(v) ((v)>>12)
#define USE_NORMAL false

std::vector<int> gather_samples(int x, int y, C2C *c2c) {
    // Grab samples and combine them into a list
    auto xy = std::make_pair(x, y);
    xy.first--;
    xy.second--;
    auto c1 = c2c->data[xy];
    xy.first++;
    auto c2 = c2c->data[xy];
    xy.first++;
    auto c3 = c2c->data[xy];
    xy.second++;
    auto c4 = c2c->data[xy];
    xy.second++;
    auto c5 = c2c->data[xy];
    xy.first--;
    auto c6 = c2c->data[xy];
    xy.first--;
    auto c7 = c2c->data[xy];
    xy.second--;
    auto c8 = c2c->data[xy];

    auto samples = std::vector<int>();
    samples.reserve(c1.size() + c2.size() + c3.size() + c4.size() +
                    c5.size() + c6.size() + c7.size() + c8.size());
    samples.insert(samples.end(), c1.begin(), c1.end());
    samples.insert(samples.end(), c2.begin(), c2.end());
    samples.insert(samples.end(), c3.begin(), c3.end());
    samples.insert(samples.end(), c4.begin(), c4.end());
    samples.insert(samples.end(), c5.begin(), c5.end());
    samples.insert(samples.end(), c6.begin(), c6.end());
    samples.insert(samples.end(), c7.begin(), c7.end());
    samples.insert(samples.end(), c8.begin(), c8.end());

    return samples;
}

std::vector<int> gather_samples_bit(int x, int y, BITMAP *a) {
    std::vector<int> samples;
    for (int dy = 0; dy < patch_w; dy++) {
        int *arow = &(*a)[y + dy][x];
        for (int dx = 0; dx < patch_w; dx++) {
            int ac = arow[dx];
            samples.push_back(ac);
        }
    }
    return samples;
}

Histogram create_hist(int x, int y, BITMAP *a) {
    return Histogram(gather_samples_bit(x, y, a));
}

Histogram create_hist(int x, int y, C2C *c2c) {
    return Histogram(gather_samples(x, y, c2c));
}

float scalar_product(const Histogram &hist_a, const Histogram &hist_b) {
    float product = 0.0f;
    for (unsigned int i = 0; i < hist_a.size(); i++) {
        product += hist_a[i] * hist_b[i];
    }
    return (1.0f - (0.5f * (product + 1.0f)));
}

int dist_scalar(BITMAP *a, C2C *data, int ax, int ay, int bx, int by, int cutoff = INT_MAX) {
    // Sample in region and create histogram
    // Create a vector with n new points in the neighborhood of ax/ay to compare
    // Default is 10 right now but this can become configurable
    //std::random_device rand_dev;
    //std::mt19937 generator(rand_dev());
    //std::uniform_int_distribution<int> xdistr(x - patch_w, x + patch_w);
    //std::uniform_int_distribution<int> ydistr(y - patch_w, y + patch_w);

    // TODO Need to make this a distribution from the BITMAP
    auto histogram_a = create_hist(ax, ay, a);
    auto histogram_b = create_hist(bx, by, data);

    // Calculate scalar product
    auto scalar = scalar_product(histogram_a, histogram_b);
    int ans = static_cast<int>(scalar * INT_MAX);
    //std::cout << "Ans: " << ans << std::endl;

    //float min_dist = std::numeric_limits<float>::max();
    //std::pair<int, int> min_xy;
    //unsigned int num_samples = 10;
    //for (unsigned int i = 0; i < num_samples; i++) {
    //  int rand_x = xdistr(generator);
    //  int rand_y = ydistr(generator);
    //  // TODO Remove these hard-coded values
    //  rand_x = std::clamp(rand_x, 0, 800);
    //  rand_y = std::clamp(rand_y, 0, 800);
    //  // Checking this since we don't want to compare
    //  // against the point being searched for
    //  if (rand_x != x && rand_y != y) {
    //    auto histogram_b = create_hist(rand_x, rand_y, data);
    //    auto scalar_prod = scalar_product(histogram_a, histogram_b);
    //    if (scalar_prod < min_dist) {
    //      min_dist = scalar_prod;
    //      min_xy = std::make_pair(rand_x, rand_y);
    //    }
    //  }
    //}

    // Calculate RGB value for candidate pixel out of dump
    //auto samples_a = gather_samples_bit(ax, ay, a);
    //int r_a = 0;
    //int g_a = 0;
    //int b_a = 0;
    //for (auto sample : samples_a) {
    //  r_a += sample & 255;
    //  g_a += (sample >> 8) & 255;
    //  b_a += sample >> 16;
    //}
    //r_a = static_cast<int>(r_a / samples_a.size());
    //g_a = static_cast<int>(g_a / samples_a.size());
    //b_a = static_cast<int>(b_a / samples_a.size());

    //auto samples = gather_samples(bx, by, data);
    //int r_b = 0;
    //int g_b = 0;
    //int b_b = 0;
    //for (auto sample : samples) {
    //  r_b += sample & 255;
    //  g_b += (sample >> 8) & 255;
    //  b_b += sample >> 16;
    //}
    //r_b = static_cast<int>(r_b / samples.size());
    //g_b = static_cast<int>(g_b / samples.size());
    //b_b = static_cast<int>(b_b / samples.size());

    //int dr = (r_a & 255) - (r_b & 255);
    //int dg = ((r_a >> 8) & 255) - ((r_b >> 8) & 255);
    //int db = (r_a >> 16) - (r_b >> 16);
    //int ans = (dr * dr) + (dg * dg) + (db * db);

    if (ans >= cutoff) {
        return cutoff;
    }

    return ans;
}

/* Measure distance between 2 patches with upper left corners (ax, ay) and (bx, by), terminating early if we exceed a cutoff distance.
   You could implement your own descriptor here. */
int dist(BITMAP *a, BITMAP *b, int ax, int ay, int bx, int by, int cutoff = INT_MAX) {
    int ans = 0;
    for (int dy = 0; dy < patch_w; dy++) {
        int *arow = &(*a)[ay + dy][ax];
        int *brow = &(*b)[by + dy][bx];
        for (int dx = 0; dx < patch_w; dx++) {
            int ac = arow[dx];
            int bc = brow[dx];
            int dr = (ac & 255) - (bc & 255);
            int dg = ((ac >> 8) & 255) - ((bc >> 8) & 255);
            int db = (ac >> 16) - (bc >> 16);
            ans += dr * dr + dg * dg + db * db;
        }
        if (ans >= cutoff) {
            return cutoff;
        }
    }
    return ans;
}

void
improve_guess(BITMAP *a, BITMAP *b, C2C *data, int ax, int ay, int &xbest, int &ybest, int &dbest, int bx, int by) {
    int d = 0;
    if (USE_NORMAL) {
        d = dist(a, b, ax, ay, bx, by, dbest);
    } else {
        d = dist_scalar(a, data, ax, by, bx, by, dbest);
    }
    if (d < dbest) {
        dbest = d;
        xbest = bx;
        ybest = by;
    }
}

/* Match image a to image b, returning the nearest neighbor field mapping a => b coords, stored in an RGB 24-bit image as (by<<12)|bx. */
void patchmatch(BITMAP *a, BITMAP *b, C2C *data, BITMAP *&ann, BITMAP *&annd) {
    /* Initialize with random nearest neighbor field (NNF). */
    ann = new BITMAP(a->w, a->h);
    annd = new BITMAP(a->w, a->h);
    int aew = a->w - patch_w + 1, aeh =
            a->h - patch_w + 1;       /* Effective width and height (possible upper left corners of patches). */
    int bew = b->w - patch_w + 1, beh = b->h - patch_w + 1;
    memset(ann->data, 0, sizeof(int) * a->w * a->h);
    memset(annd->data, 0, sizeof(int) * a->w * a->h);
    for (int ay = 0; ay < aeh; ay++) {
        for (int ax = 0; ax < aew; ax++) {
            int bx = rand() % bew;
            int by = rand() % beh;
            (*ann)[ay][ax] = XY_TO_INT(bx, by);
            if (USE_NORMAL) {
                (*annd)[ay][ax] = dist(a, b, ax, ay, bx, by);
            } else {
                (*annd)[ay][ax] = dist_scalar(a, data, ax, ay, bx, by);
            }
        }
    }

    std::cout << "Init done" << std::endl;

    for (int iter = 0; iter < pm_iters; iter++) {
        /* In each iteration, improve the NNF, by looping in scanline or reverse-scanline order. */
        int ystart = 0, yend = aeh, ychange = 1;
        int xstart = 0, xend = aew, xchange = 1;
        if (iter % 2 == 1) {
            xstart = xend - 1;
            xend = -1;
            xchange = -1;
            ystart = yend - 1;
            yend = -1;
            ychange = -1;
        }
        for (int ay = ystart; ay != yend; ay += ychange) {
            for (int ax = xstart; ax != xend; ax += xchange) {
                /* Current (best) guess. */
                int v = (*ann)[ay][ax];
                int xbest = INT_TO_X(v), ybest = INT_TO_Y(v);
                int dbest = (*annd)[ay][ax];

                /* Propagation: Improve current guess by trying instead correspondences from left and above (below and right on odd iterations). */
                if ((unsigned) (ax - xchange) < (unsigned) aew) {
                    int vp = (*ann)[ay][ax - xchange];
                    int xp = INT_TO_X(vp) + xchange, yp = INT_TO_Y(vp);
                    if ((unsigned) xp < (unsigned) bew) {
                        improve_guess(a, b, data, ax, ay, xbest, ybest, dbest, xp, yp);
                    }
                }

                if ((unsigned) (ay - ychange) < (unsigned) aeh) {
                    int vp = (*ann)[ay - ychange][ax];
                    int xp = INT_TO_X(vp), yp = INT_TO_Y(vp) + ychange;
                    if ((unsigned) yp < (unsigned) beh) {
                        improve_guess(a, b, data, ax, ay, xbest, ybest, dbest, xp, yp);
                    }
                }

                /* Random search: Improve current guess by searching in boxes of exponentially decreasing size around the current best guess. */
                int rs_start = rs_max;
                if (rs_start > MAX(b->w, b->h)) { rs_start = MAX(b->w, b->h); }
                for (int mag = rs_start; mag >= 1; mag /= 2) {
                    /* Sampling window */
                    int xmin = MAX(xbest - mag, 0), xmax = MIN(xbest + mag + 1, bew);
                    int ymin = MAX(ybest - mag, 0), ymax = MIN(ybest + mag + 1, beh);
                    int xp = xmin + rand() % (xmax - xmin);
                    int yp = ymin + rand() % (ymax - ymin);
                    improve_guess(a, b, data, ax, ay, xbest, ybest, dbest, xp, yp);
                }

                (*ann)[ay][ax] = XY_TO_INT(xbest, ybest);
                (*annd)[ay][ax] = dbest;
            }
        }
    }
}

int main(int argc, char *argv[]) {
    argc--;
    argv++;
    if (argc != 5) {
        fprintf(stderr, "pm_minimal a b ann annd\n"
                        "Given input images a, b outputs nearest neighbor field 'ann' mapping a => b coords, and the squared L2 distance 'annd'\n"
                        "These are stored as RGB 24-bit images, with a 24-bit int at every pixel. For the NNF we store (by<<12)|bx.");
        exit(1);
    }
    printf("Loading input images\n");
    BITMAP *a = load_bitmap(argv[0]);
    BITMAP *b = load_bitmap(argv[1]);
    C2C *data = load_dump(argv[2]);
    BITMAP *ann = NULL, *annd = NULL;
    printf("Running PatchMatch\n");
    patchmatch(a, b, data, ann, annd);
    printf("Saving output images\n");
    save_bitmap(ann, argv[3]);
    save_bitmap(annd, argv[4]);
    return 0;
}
