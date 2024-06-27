static class main{
		static void Main(){
			int i = 1;
			while(i+1>i){i++;}
			System.Console.Write($"My max int = {i}\n");
			i = 1;
			while(i-1<i){i--;}
			System.Console.Write($"My min int = {i}\n");
			double doubleMEpsilon=1;while(doubleMEpsilon+1!=1){doubleMEpsilon/=2;} doubleMEpsilon*=2;
			float floatMEpsilon=1;while(floatMEpsilon+1!=1){floatMEpsilon/=2;} floatMEpsilon*=2;
			System.Console.Write($"Double machine epsilon = {doubleMEpsilon}\n");
			System.Console.Write($"Float machine epsilon = {floatMEpsilon}\n");
			System.Console.Write($"Theoretical Double machine epsilon = {System.Math.Pow(2, -52)}\n");
			System.Console.Write($"Theoretical float machine epsilon = {System.Math.Pow(2, -23)}\n");
			double tiny = doubleMEpsilon/2;
			double lessTiny1 = 1+tiny+tiny;
			double lessTiny2 = tiny+tiny+1;
			System.Console.Write($"1+tiny+tiny=tiny+tiny+1 ? {lessTiny1==lessTiny2}\n");
			System.Console.Write($"1+tiny+tiny>=tiny+tiny+1 ? {lessTiny1>lessTiny2}\n");
			System.Console.Write($"1+tiny+tiny<tiny+tiny+1 ? {lessTiny1<lessTiny2}\n");
			System.Console.Write("\n Test of double comparisons\n");
			double d1 = 0.1 + 0.1 + 0.1 + 0.1 + 0.1 + 0.1 + 0.1 + 0.1;
			double d2 = 8*0.1;
			System.Console.Write($"d1 = {d1:e15}\n");
			System.Console.Write($"d2 = {d2:e15}\n");
			System.Console.Write($"d1 == d2 ? => {d1==d2}\n");
			System.Console.Write($"d1 ~ d2 ? => {approx(d1, d2)}\n");
		}
		public static bool approx(double a, double b, double acc = 1e-9, double eps=1e-9)
		{
			bool condition1 = System.Math.Abs(a - b) < acc;
			bool condition2 = System.Math.Abs(a - b)/System.Math.Max(System.Math.Abs(a), System.Math.Abs(b)) < eps;
			return condition1 | condition2;
		}
}
