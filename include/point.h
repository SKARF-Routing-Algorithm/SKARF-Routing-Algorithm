#pragma once
struct Point {
   private:
    double lat, lon;

   public:
    Point(double _lon, double _lat) {
        lat = _lat;
        lon = _lon;
    };
    double get_lat() { return lat; }
    double get_lon() { return lon; }
    void set_lat(double _lat) { lat = _lat; };
    void set_lon(double _lon) { lon = _lon; };

    bool operator<(const Point& q) const {
        if (lat == q.lat) return lon < q.lon;
        return lat < q.lat;
    }
};