//FILE: SanityCheck.cs

using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;


public static class SanityCheck
{
    public static void AssertFailed()
    {
        throw new Exception("User-defined exception");
    }
}
