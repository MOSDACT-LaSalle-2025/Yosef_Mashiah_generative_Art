/* 
Name: Yosef Mashiah
Date: 31.8.2025
Description: 
A Processing sketch that draws a grid of square cells. 
Each cell displays a gradient from white to a unique color, with a rotating hand in the center. 
When the hand of one cell touches the hand of a neighboring cell, 
that cell reverses its rotation direction and changes to a new color. 
Place of production: Palastine
*/



int cols = 10;
int rows = 10;
float cellSize;

float[][] angles;
float[][] omega;
int[][] cooldown;
color[][] baseColors;

float[][] cx, cy;
float[][] ex, ey;

float r;
float TOUCH_EPS;

void setup() {
  size(600, 600, P2D);
  cellSize = floor((float)width / cols);

  angles     = new float[cols][rows];
  omega      = new float[cols][rows];
  cooldown   = new int[cols][rows];
  baseColors = new color[cols][rows];

  cx = new float[cols][rows];
  cy = new float[cols][rows];
  ex = new float[cols][rows];
  ey = new float[cols][rows];

  r = cellSize * 0.5 - 1;
  TOUCH_EPS = max(2.0, cellSize * 0.1);

  for (int i = 0; i < cols; i++) {
    for (int j = 0; j < rows; j++) {
      angles[i][j] = random(TWO_PI);
      float speed = random(0.01, 0.03) * (random(1) < 0.5 ? -1 : 1);
      omega[i][j] = speed;
      cooldown[i][j] = 0;
      baseColors[i][j] = color(random(255), random(255), random(255));

      cx[i][j] = i * cellSize + cellSize * 0.5;
      cy[i][j] = j * cellSize + cellSize * 0.5;
    }
  }

  noStroke();
  smooth(8);
}

void draw() {
  background(255); // נצבע הכל בלבן כל פריים

  for (int i = 0; i < cols; i++) {
    for (int j = 0; j < rows; j++) {
      // ציור גרדיאנט בריבוע התא (לבן→צבע ייחודי)
      drawSquareGradient(i * cellSize, j * cellSize, cellSize, baseColors[i][j]);

      // ציור המחוג – מסובב סביב מרכז התא
      pushMatrix();
      translate(cx[i][j], cy[i][j]);
      rotate(angles[i][j]);
      stroke(0);
      strokeWeight(2);
      line(0, 0, r, 0);
      popMatrix();

      ex[i][j] = cx[i][j] + r * cos(angles[i][j]);
      ey[i][j] = cy[i][j] + r * sin(angles[i][j]);
    }
  }

  // בדיקת נגיעות
  for (int i = 0; i < cols; i++) {
    for (int j = 0; j < rows; j++) {
      if (i + 1 < cols) {
        if (segmentsClose(cx[i][j], cy[i][j], ex[i][j], ey[i][j],
                          cx[i+1][j], cy[i+1][j], ex[i+1][j], ey[i+1][j], TOUCH_EPS)) {
          flipIfReady(i, j);
          flipIfReady(i+1, j);
        }
      }
      if (j + 1 < rows) {
        if (segmentsClose(cx[i][j], cy[i][j], ex[i][j], ey[i][j],
                          cx[i][j+1], cy[i][j+1], ex[i][j+1], ey[i][j+1], TOUCH_EPS)) {
          flipIfReady(i, j);
          flipIfReady(i, j+1);
        }
      }
    }
  }

  for (int i = 0; i < cols; i++) {
    for (int j = 0; j < rows; j++) {
      angles[i][j] += omega[i][j];
      if (cooldown[i][j] > 0) cooldown[i][j]--;
    }
  }
}

// גרדיאנט מרובע מלבן לצבע
void drawSquareGradient(float x, float y, float size, color base) {
  noStroke();
  for (int yy = 0; yy < int(size); yy++) {
    float t = map(yy, 0, size - 1, 0, 1);
    color c = lerpColor(color(255), base, t);
    fill(c);
    rect(x, y + yy, size, 1.1);
  }
}

// הפיכת כיוון + צבע
void flipIfReady(int i, int j) {
  if (cooldown[i][j] <= 0) {
    omega[i][j] *= -1;
    baseColors[i][j] = color(random(255), random(255), random(255));
    cooldown[i][j] = 10;
  }
}

boolean segmentsClose(float ax, float ay, float bx, float by,
                      float cx_, float cy_, float dx, float dy, float eps) {
  return segmentDistance(ax, ay, bx, by, cx_, cy_, dx, dy) <= eps;
}

float segmentDistance(float ax, float ay, float bx, float by,
                      float cx, float cy, float dx, float dy) {
  if (ax == bx && ay == by) return pointToSegment(cx, cy, dx, dy, ax, ay);
  if (cx == dx && cy == dy) return pointToSegment(ax, ay, bx, by, cx, cy);
  float ux = bx - ax, uy = by - ay;
  float vx = dx - cx, vy = dy - cy;
  float wx = ax - cx, wy = ay - cy;
  float a = ux*ux + uy*uy;
  float b = ux*vx + uy*vy;
  float c = vx*vx + vy*vy;
  float d = ux*wx + uy*wy;
  float e = vx*wx + vy*wy;
  float D = a*c - b*b;
  float sc, sN, sD = D;
  float tc, tN, tD = D;
  if (D < 1e-8) { sN = 0; sD = 1; tN = e; tD = c; }
  else {
    sN = (b*e - c*d);
    tN = (a*e - b*d);
    if (sN < 0) { sN = 0; tN = e; tD = c; }
    else if (sN > sD) { sN = sD; tN = e + b; tD = c; }
  }
  if (tN < 0) {
    tN = 0;
    if (-d < 0) sN = 0;
    else if (-d > a) sN = sD;
    else { sN = -d; sD = a; }
  } else if (tN > tD) {
    tN = tD;
    if ((-d + b) < 0) sN = 0;
    else if ((-d + b) > a) sN = sD;
    else { sN = (-d + b); sD = a; }
  }
  sc = (abs(sN) < 1e-8 ? 0 : sN / sD);
  tc = (abs(tN) < 1e-8 ? 0 : tN / tD);
  float dx_ = wx + sc * ux - tc * vx;
  float dy_ = wy + sc * uy - tc * vy;
  return sqrt(dx_*dx_ + dy_*dy_);
}

float pointToSegment(float ax, float ay, float bx, float by, float px, float py) {
  float vx = bx - ax, vy = by - ay;
  float wx = px - ax, wy = py - ay;
  float vLen2 = vx*vx + vy*vy;
  float t = (vLen2 < 1e-8) ? 0 : (wx*vx + wy*vy) / vLen2;
  t = constrain(t, 0, 1);
  float qx = ax + t * vx;
  float qy = ay + t * vy;
  float dx = px - qx, dy = py - qy;
  return sqrt(dx*dx + dy*dy);
}
