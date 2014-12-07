// FILE: LocalExtensions.cs

using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.IO;
using System.Globalization;

public static class LocalExtensions
{

	public static Random RNG = new Random();

    public static void WriteSep(this StreamWriter sw, string sep, bool trailingSep, params object[] objs)
    {
        for (int i = 0; i < objs.Length; i++)
        {
            sw.Write(objs[i].ToString());
            if (trailingSep || i != objs.Length - 1) sw.Write(sep);
        }
    }

    public static void WriteSepLine(this StreamWriter sw, string sep, params object[] objs)
    {
        WriteSep(sw, sep, false, objs);
        sw.WriteLine();
    }

    public static string ToPercentage(this double d)
    {
        return (d * 100).ToString("N1") + "%";
    }
}