//FILE: Maths.cs

using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Drawing.Drawing2D;
using System.Drawing;

public static class Maths
{

    public static double ToRadians(this double val)
    {
        return (Math.PI / 180) * val;
    }

    public static double ToDegrees(this double val)
    {
        return (180 * val) / Math.PI;
    }

    // Mercator projection -- from http://stackoverflow.com/questions/18838915/convert-lat-lon-to-pixel-coordinate
    public static PointF ProjectToMercator(Coord latlonInDegrees, MercatorInfo mi)
    {
        Coord latLonInRads = new Coord(latlonInDegrees.Lat.ToRadians(), latlonInDegrees.Lon.ToRadians());
        double x = latLonInRads.Lon;
        double y = MercatorInfo.ToMercatorY(latLonInRads.Lat);
        x = (x - mi.west) * mi.xFactor;
        y = (mi.ymax - y) * mi.yFactor; // y points south
        return new PointF((float)x, (float)y);
    }

    public static Coord ProjectBackFromMercator(PointF pt, MercatorInfo mi)
    {
        double lon = (pt.X / mi.xFactor) + mi.west;
        double lat = 2 * Math.Atan(Math.Exp(mi.ymax - (pt.Y / mi.yFactor))) - .5 * Math.PI;

        return new Coord(lat.ToDegrees(), lon.ToDegrees());
    }

    public class MercatorInfo
    {
        public double ymin;
        public double ymax;
        public double xFactor;
        public double yFactor;
        public double west;

        public MercatorInfo(double bottomLatInRads, double topLatInRads, double leftLonInRads, double rightLonInRads, double aspectWidth, double aspectHeight)
        {
            ymin = ToMercatorY(bottomLatInRads);
            ymax = ToMercatorY(topLatInRads);
            xFactor = aspectWidth / (rightLonInRads - leftLonInRads);
            yFactor = aspectHeight / (ymax - ymin);
            west = leftLonInRads;
        }

        public static double ToMercatorY(double latInRads)
        {
            return Math.Log(Math.Tan((latInRads * .5) + (Math.PI * .25)));
        }
    }

    #region Distance measures

    internal static double HaversineDistance(double lat1, double lon1, double lat2, double lon2)
    {
        double R = 6371; // km
        double dLat = (lat2 - lat1) * (Math.PI / 180);
        double dLon = (lon2 - lon1) * (Math.PI / 180);
        lat1 = lat1 * (Math.PI / 180);
        lat2 = lat2 * (Math.PI / 180);

        double a = Math.Sin(dLat / 2) * Math.Sin(dLat / 2) +
                Math.Sin(dLon / 2) * Math.Sin(dLon / 2) * Math.Cos(lat1) * Math.Cos(lat2);
        double c = 2 * Math.Atan2(Math.Sqrt(a), Math.Sqrt(1 - a));
        return R * c;
    }

    internal static double GISDistance(double lat1, double lon1, double lat2, double lon2)
    {
        const double R = 6371; // km
        lat1 = lat1 * (Math.PI / 180);
        lat2 = lat2 * (Math.PI / 180);
        lon1 = lon1 * (Math.PI / 180);
        lon2 = lon2 * (Math.PI / 180);

        return Math.Acos(Math.Sin(lat1) * Math.Sin(lat2) +
                          Math.Cos(lat1) * Math.Cos(lat2) *
                          Math.Cos(lon2 - lon1)) * R;
    }

    internal static double EucDistance(double x1, double y1, double x2, double y2)
    {
        double xDiff = (x2 - x1);
        double yDiff = (y2 - y1);
        double xDiffSquared = (xDiff * xDiff);
        double yDiffSquared = (yDiff * yDiff);
        if (xDiffSquared < 0 || yDiffSquared < 0) SanityCheck.AssertFailed();
        return Math.Sqrt(xDiffSquared + yDiffSquared);
    }


    #endregion

}
