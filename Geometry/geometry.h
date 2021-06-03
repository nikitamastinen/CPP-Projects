//
// Created by Nikita Mastinen on 15.11.2020.
//

#include <utility>
#include <iostream>
#include <vector>
#include <cmath>

bool isDoublesEqual(const double a, const double b, double precision = 1e-6) {
  return std::abs(a - b) <= precision;
}

struct Point;
struct Vector;

struct Point {
  double x = 0.0, y = 0.0;

  Point() = default;

  explicit Point(const Vector& a);

  Point(const double x, const double y): x(x), y(y) {};

  bool operator==(const Point& other) const {
    return isDoublesEqual(x, other.x) && isDoublesEqual(y, other.y);
  }

  bool operator!=(const Point& other) const {
    return !isDoublesEqual(x, other.x) || !isDoublesEqual(y, other.y);
  }
};

struct Vector {
  double x = 0.0, y = 0.0;

  Vector() = default;

  Vector(const double x, const double y): x(x), y(y) {};

  Vector(const Point& a, const Point& b): x(b.x - a.x), y(b.y - a.y) {};

  explicit Vector(const Point& a): x(a.x), y(a.y) {};

  Vector& normalize() {
    double k = std::pow(x * x + y * y, 0.5);
    x /= k, y /= k;
    return *this;
  }

  double length() const {
    return std::pow(x * x + y * y, 0.5);
  }

  bool operator==(const Vector& other) const {
    return isDoublesEqual(x, other.x) && isDoublesEqual(y, other.y);
  }

  Vector rotate(const double phi) {
    double x_copy = x;
    x = x * cos(phi) - y * sin(phi);
    y = x_copy * sin(phi) + y * cos(phi);
    return *this;
  }

  bool operator!=(const Vector& other) const {
    return !isDoublesEqual(x, other.x) || !isDoublesEqual(y, other.y);
  }
};

Vector operator*(const Vector& a, const double k) {
  return {a.x * k, a.y * k};
}

double multiply(const Vector& a, const Vector& b) {
  return a.x * b.x + a.y * b.y;
}

Vector operator+(const Vector& a, const Vector& b) {
  return {a.x + b.x, a.y + b.y};
}

Vector operator-(const Vector& a, const Vector& b) {
  return {a.x - b.x, a.y - b.y};
}

double operator*(const Vector& a, const Vector& b) {
  return a.x * b.y - b.x * a.y;
}

Vector projection(const Vector& a, const Vector& b) {
  return b * (multiply(a, b) / b.length() / b.length());
}

Point::Point(const Vector& a): x(a.x), y(a.y) {}

class Line {
public:
  double a = 1.0, b = 0.0, c = 0.0;

  Line() = default;

  Line(const double a, const double b, const double c): a(a), b(b), c(c) {};

  Line(const Point& u, const Point& v): a(v.y - u.y),
                                        b(u.x - v.x),
                                        c(-u.x * v.y + v.x * u.y) {};

  explicit Line(const Point& u, const Vector& t): Line(u, Point(static_cast<Vector>(u) + t)) {};

  Line(const double k, const double b): a(k), b(-1.0), c(b) {};

  Line(const Point& a, const double & k): Line(k, a.y - a.x * k) {};
};

bool operator==(const Line& l, const Line& m) {
  double det_ab = l.a * m.b - m.a * l.b;
  double det_ac = l.a * m.c - m.a * l.c;
  return isDoublesEqual(det_ab, 0.0) && isDoublesEqual(det_ac, 0.0);
}

Point intersect(Line l, Line m) {
  double det = l.a * m.b - m.a * l.b;
  return {(m.c * l.b - l.c * m.b) / det, (-m.c * l.a + l.c * m.a) / det};
}

class Shape {
public:
  virtual ~Shape() = default;

  virtual double perimeter() const = 0;

  virtual double area() const = 0;

  virtual bool operator==(const Shape& other) const = 0;

  virtual bool operator!=(const Shape& other) const = 0;

  virtual bool isCongruentTo(const Shape& other) const = 0;

  virtual bool isSimilarTo(const Shape& other) const = 0;

  virtual bool containsPoint(const Point& point) const = 0;

  virtual void rotate(const Point& center, double angle) = 0;

  virtual void reflex(const Point& center) = 0;

  virtual void reflex(const Line& axis) = 0;

  virtual void scale(Point center, double coefficient) = 0;
};

class Polygon: public Shape {
protected:
  std::vector<Point> vertices;
public:
  Polygon() = default;

  explicit Polygon(const std::vector<Point>& vertices): vertices(vertices) {};

  Polygon(const std::initializer_list<Point>& l) : vertices(l) {};

  size_t verticesCount() const {
    return vertices.size();
  }

  const std::vector<Point>& getVertices() const {
    return vertices;
  }

  bool isConvex() const {
    bool more_zero = false, less_zero = false;
    for (size_t i = 0; i < vertices.size(); i++) {
      if (Vector(vertices[i], vertices[(i + 1) % vertices.size()]) *
          Vector(vertices[(i + 1) % vertices.size()],
                 vertices[(i + 2) % vertices.size()]) > 0) {
        more_zero = true;
      } else {
        less_zero = true;
      }
    }
    return !more_zero || !less_zero;
  }

  double perimeter() const override {
    double result = Vector(vertices.back(), vertices.front()).length();
    for (size_t i = 0; i + 1 < vertices.size(); ++i) {
      result += Vector(vertices[i], vertices[i + 1]).length();
    }
    return false;
  }

  double area() const override {
    double res = 0;
    for (size_t i = 1; i + 1 < vertices.size(); i++) {
      res += Vector(vertices[0], vertices[i]) * Vector(vertices[0], vertices[i + 1]);
    }
    return std::abs(res / 2);
  }

  bool operator==(const Shape& other) const override {
    if (dynamic_cast<const Polygon*>(&other) == nullptr) {
      return false;
    }
    const auto& another = dynamic_cast<const Polygon&>(other);
    if (another.vertices.size() != vertices.size()) return false;
    for (size_t i = 0; i < 2 * vertices.size(); ++i) {
      for (size_t j = 0; j < vertices.size(); ++j) {
        if (vertices[(j + i) % vertices.size()] != another.vertices[j]) break;
        if (j + 1 == vertices.size()) return true;
      }
      for (size_t j = vertices.size() - 1; ; --j) {
        if (vertices[((int)vertices.size() - j - 1 + i) %
                     vertices.size()] != another.vertices[j]) break;
        if (j == 0) return true;
      }
    }
    return false;
  }

  bool operator!=(const Shape& other) const override {
    return !(*this == other);
  }

  bool isCongruentTo(const Shape& other) const override {
    if (dynamic_cast<const Polygon*>(&other) == nullptr) {
      return false;
    }
    const auto& another = dynamic_cast<const Polygon&>(other);
    if (another.vertices.size() != vertices.size()) return false;
    Point this_centroid, other_centroid;
    for (size_t i = 0; i < verticesCount(); i++) {
      this_centroid = static_cast<Point>(static_cast<Vector>(this_centroid) +
                                         static_cast<Vector>(vertices[i]));
      other_centroid = static_cast<Point>(static_cast<Vector>(other_centroid) +
                                          static_cast<Vector>(another.vertices[i]));
    }
    this_centroid = static_cast<Point>(static_cast<Vector>(this_centroid) *
                                       (1.0 / vertices.size()));
    other_centroid = static_cast<Point>(static_cast<Vector>(other_centroid) *
                                        (1.0 / another.vertices.size()));
    Polygon this_copy, other_copy;
    for (size_t i = 0; i < vertices.size(); i++) {
      this_copy.vertices.emplace_back(Vector(static_cast<Vector>(vertices[i]) -
                                             static_cast<Vector>(this_centroid)).length(),
                                      std::abs(Vector(vertices[i], vertices[(i + 1) % vertices.size()]) *
                                               Vector(vertices[i], vertices[(i - 1 + vertices.size()) % vertices.size()])));
      other_copy.vertices.emplace_back(Vector(static_cast<Vector>(another.vertices[i]) -
                                              static_cast<Vector>(other_centroid)).length(),
                                       std::abs(Vector(another.vertices[i], another.vertices[(i + 1) % vertices.size()]) *
                                                Vector(another.vertices[i],
                                                       another.vertices[(i - 1 + vertices.size()) % vertices.size()])));
    }
    return this_copy == other_copy;
  }

  bool isSimilarTo(const Shape& other) const override {
    if (dynamic_cast<const Polygon*>(&other) == nullptr) {
      return false;
    }
    const auto& another = dynamic_cast<const Polygon&>(other);
    if (another.vertices.size() != vertices.size()) return false;
    double this_minimum = 1e100, other_minimum = 1e100;
    for (size_t i = 0; i < vertices.size(); i++) {
      this_minimum = std::min(this_minimum,
                              Vector(vertices[i], vertices[(i + 1) % vertices.size()]).length());
      other_minimum = std::min(other_minimum,
                               Vector(another.vertices[i], another.vertices[(i + 1) % vertices.size()]).length());
    }
    Polygon copy = another;
    copy.scale(Point(), this_minimum / other_minimum);
    return copy.isCongruentTo(*this);
  }

  bool containsPoint(const Point& point) const override {
    Line shot = Line(point, Point(11, 111));
    int result = 0;
    for (size_t i = 0; i < vertices.size(); i++) {
      Point intersection = intersect(shot, Line(vertices[i],
                                                vertices[(i + 1) % vertices.size()]));
      if ((static_cast<Vector>(intersection) - static_cast<Vector>(point)).x >= 0 &&
          multiply(Vector(intersection, vertices[i]),
                   Vector(intersection, vertices[(i + 1) % vertices.size()])) < 0) {
        result++;
      }
    }
    return result % 2 == 1;
  }

  void rotate(const Point& center, double angle) override {
    for (auto& vertice : vertices) {
      vertice = static_cast<Point>(static_cast<Vector>(center) +
                                   Vector(static_cast<Vector>(vertice) -
                                          static_cast<Vector>(center)).rotate(angle));
    }
  }

  void reflex(const Point& center) override {
    for (auto& vertice : vertices) {
      vertice = static_cast<Point>(static_cast<Vector>(center) +
                                   Vector(static_cast<Vector>(vertice) -
                                          static_cast<Vector>(center)) * -1.0);
    }
  }

  void reflex(const Line& axis) override {
    Vector direction = Vector(-axis.b, axis.a);
    Point any = Point(axis.c * axis.a / direction.length(),
                      axis.c * axis.b / direction.length());
    for (auto& vertice : vertices) {
      Point prev_vertice = vertice;
      vertice = static_cast<Point>(static_cast<Vector>(any) +
                                   projection(Vector(any, vertice), direction));
      vertice = static_cast<Point>(static_cast<Vector>(vertice) +
                                   Vector(prev_vertice, vertice));
    }
  }

  void scale(Point center, double coefficient) override {
    for (auto& vertice : vertices) {
      vertice = static_cast<Point>(static_cast<Vector>(center) +
                                   Vector(static_cast<Vector>(vertice) -
                                          static_cast<Vector>(center)) * coefficient);
    }
  }
};

class Ellipse: public Shape {
protected:
  std::pair<Point, Point> focus;
  double sum_of_distances = 0.0;
public:
  Ellipse() = default;

  Ellipse(const Point& f1, const Point& f2, const double sum_of_distances):
      focus({f1, f2}), sum_of_distances(sum_of_distances) {};

  const std::pair<Point, Point>& focuses() const {
    return focus;
  }

  std::pair<Line, Line> directrices() const {
    Vector direction = Vector(center(), focus.first).normalize();
    Vector normal = direction.rotate(acos(-1) / 2);
    double coefficient = sum_of_distances / 2 * eccentricity();
    return {Line(static_cast<Point>(static_cast<Vector>(center()) +
                                    static_cast<Vector>(direction) * coefficient), normal),
            Line(static_cast<Point>(static_cast<Vector>(center()) +
                                    static_cast<Vector>(direction * -coefficient)), normal)};
  }

  double eccentricity() const {
    return std::abs(Vector(static_cast<Vector>(focus.first) -
                           static_cast<Vector>(center())).length() / sum_of_distances * 2);
  }

  Point center() const {
    return static_cast<Point>((static_cast<Vector>(focus.first) +
                               static_cast<Vector>(focus.second)) * 0.5);
  }

  double perimeter() const override {
    double a = sum_of_distances / 2;
    double b = pow(a * a - eccentricity() * eccentricity() * a * a, 0.5);
    return acos(-1) * (3 * (a + b) - std::pow((3 * a + b) *  (a + 3 * b), 0.5));
  }

  double area() const override {
    double a = sum_of_distances / 2;
    return acos(-1) * a *
           std::pow(a * a - eccentricity() * eccentricity() * a * a, 0.5);
  }

  bool operator==(const Shape& other) const override {
    if (dynamic_cast<const Ellipse*>(&other) == nullptr) {
      return false;
    }
    const auto& another = dynamic_cast<const Ellipse&>(other);
    return isDoublesEqual(sum_of_distances, another.sum_of_distances) &&
           (focus == another.focus || (focus.first == another.focus.second &&
                                       focus.second == another.focus.first));
  }

  bool operator!=(const Shape& other) const override {
    return !(*this == other);
  }

  bool isCongruentTo(const Shape& other) const override {
    if (dynamic_cast<const Ellipse*>(&other) == nullptr) {
      return false;
    }
    const auto& another = dynamic_cast<const Ellipse&>(other);
    return isDoublesEqual(sum_of_distances, another.sum_of_distances) &&
           isDoublesEqual(Vector(focus.first, focus.second).length(),
                          Vector(another.focus.first, another.focus.second).length());
  }

  bool isSimilarTo(const Shape& other) const override {
    if (dynamic_cast<const Ellipse*>(&other) == nullptr) {
      return false;
    }
    auto another = dynamic_cast<const Ellipse&>(other);
    double coefficient = sum_of_distances / another.sum_of_distances;
    another.scale(Point(), coefficient);
    return isCongruentTo(another);
  }

  bool containsPoint(const Point& point) const override {
    return Vector(point, focus.first).length() +
           Vector(point, focus.second).length() <= sum_of_distances;
  }

  void rotate(const Point& center, double angle) override {
    focus.first = static_cast<Point>(static_cast<Vector>(focus.first) +
                                     Vector(center, focus.first).rotate(angle));
    focus.second = static_cast<Point>(static_cast<Vector>(focus.second) +
                                      Vector(center, focus.second).rotate(angle));
  }

  void reflex(const Point& center) override {
    scale(center, -1);
  }

  void reflex(const Line& axis) override {
    Polygon polygon = {focus.first, focus.second};
    polygon.reflex(axis);
    focus = {polygon.getVertices()[0], polygon.getVertices()[1]};
  }

  void scale(Point center, double coefficient) override {
    focus.first = static_cast<Point>(static_cast<Vector>(center) +
                                     Vector(static_cast<Vector>(focus.first) - static_cast<Vector>(center)) * coefficient);
    focus.second = static_cast<Point>(static_cast<Vector>(center) +
                                      Vector(static_cast<Vector>(focus.second) - static_cast<Vector>(center)) * coefficient);
    sum_of_distances *= coefficient;
  }
};

class Circle: public Ellipse {
public:
  Circle() = default;

  Circle(Point center, double r): Ellipse(center, center, 2.0 * r) {};

  Circle(Point center, Point a): Circle(center, Vector(center, a).length()) {};

  double radius() const {
    return sum_of_distances / 2.0;
  }
};

class Rectangle: public Polygon {
public:
  Rectangle(const Point& a, const Point& c, double k) {
    if (k < 1.0) k = 1.0 / k;
    Point b = static_cast<Point>(static_cast<Vector>(a) + (Vector(a, c) *
                                                           (1.0 / std::pow(1.0 + k * k, 0.5))).
        rotate(asin(1.0 * k / std::pow(1.0 + k * k, 0.5))));
    Point d = static_cast<Point>(static_cast<Vector>(a) + (Vector(a, c) *
                                                           (k / std::pow(1.0 + k * k, 0.5))).
        rotate(asin(-1.0 / std::pow(1.0 + k * k, 0.5))));
    vertices = {a, b, c, d};
  };

  Point center() {
    return Point((vertices[0].x + vertices[2].x) / 2.0,
                 (vertices[0].y + vertices[2].y) / 2.0);
  }

  std::pair<Line, Line> diagonals() {
    return {Line(vertices[0], vertices[2]), Line(vertices[1], vertices[3])};
  }
};

class Square: public Rectangle {
public:
  Square(Point a, Point c): Rectangle(a, c, 1.0) {};

  Circle circumscribedCircle() {
    return Circle(center(),
                  std::abs(vertices[0].x - vertices[1].x) / std::pow(2.0, 0.5));
  }

  Circle inscribedCircle() {
    return Circle(center(), std::abs(vertices[0].x - vertices[1].x) / 2.0);
  }
};

class Triangle: public Polygon {
public:
  Triangle(Point a, Point b, Point c): Polygon({a, b, c}) {};

  Circle circumscribedCircle() {
    Point circumscribed_center = static_cast<Point>(static_cast<Vector>(orthocenter()) +
                                                    Vector(orthocenter(), centroid()) * 1.5);
    return Circle(circumscribed_center, vertices[0]);
  }

  Circle inscribedCircle() {
    Line first_bisector = Line(vertices[0], Point(((static_cast<Vector>(vertices[1]) -
                                                    static_cast<Vector>(vertices[0])).normalize() +
                                                   (static_cast<Vector>(vertices[2]) - static_cast<Vector>(vertices[0])).normalize()) +
                                                  static_cast<Vector>(vertices[0])));
    Line second_bisector = Line(vertices[1], Point(((static_cast<Vector>(vertices[2]) -
                                                     static_cast<Vector>(vertices[1])).normalize() +
                                                    (static_cast<Vector>(vertices[0]) - static_cast<Vector>(vertices[1])).normalize()) +
                                                   static_cast<Vector>(vertices[1])));
    Point in_center = intersect(first_bisector, second_bisector);
    Vector projector = projection(Vector(vertices[0], in_center),
                                  Vector(vertices[0], vertices[1])) + static_cast<Vector>(vertices[0]);
    return Circle(in_center, Point(projector));
  }

  Point centroid() {
    return static_cast<Point>((static_cast<Vector>(vertices[0]) +
                               static_cast<Vector>(vertices[1]) + static_cast<Vector>(vertices[2])) * (1.0 / 3.0));
  }

  Point orthocenter() {
    Line first_altitude = Line(vertices[0],
                               Point(projection(static_cast<Vector>(vertices[0]) -
                                                static_cast<Vector>(vertices[1]), static_cast<Vector>(vertices[2]) -
                                                                                  static_cast<Vector>(vertices[1])) + static_cast<Vector>(vertices[1])));
    Line second_altitude = Line(vertices[1],
                                Point(projection(static_cast<Vector>(vertices[1]) -
                                                 static_cast<Vector>(vertices[2]), static_cast<Vector>(vertices[0]) -
                                                                                   static_cast<Vector>(vertices[2])) + static_cast<Vector>(vertices[2])));
    return intersect(first_altitude, second_altitude);
  }

  Line EulerLine() {
    return Line(orthocenter(), centroid());
  }

  Circle ninePointsCircle() {
    auto result = circumscribedCircle();
    result.scale(orthocenter(), 0.5);
    return result;
  }
};