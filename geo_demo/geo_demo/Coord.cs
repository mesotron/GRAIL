// FILE: Coord.cs

using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;


/// <summary>
/// A lat/lon coordinate.
/// </summary>
public class Coord
{
    public double Lat;
    public double Lon;

    /// <summary>
    /// Constructor.
    /// </summary>
    public Coord(double lat, double lon)
    {
        this.Lat = lat;
        this.Lon = lon;
    }

    public Coord(string latlon)
    {
        string[] tokens = latlon.Split(',');
        Lat = double.Parse(tokens[0]);
        Lon = double.Parse(tokens[1]);
    }

    public override string  ToString()
    {
        return Lat.ToString() + "," + Lon.ToString();
    }

    // override object.Equals
    public override bool Equals(object obj)
    {
        //       
        // See the full list of guidelines at
        //   http://go.microsoft.com/fwlink/?LinkID=85237  
        // and also the guidance for operator== at
        //   http://go.microsoft.com/fwlink/?LinkId=85238
        //

        if (obj == null || GetType() != obj.GetType())
        {
            return false;
        }

        Coord other = (Coord)obj;
        return (this.Lat == other.Lat && this.Lon == other.Lon);
    }

    // override object.GetHashCode
    public override int GetHashCode()
    {
        return (this.Lon.ToString() + this.Lon.ToString()).GetHashCode();
    }
}
