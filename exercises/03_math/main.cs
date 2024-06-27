using System;

class main{
	static int x = 2;
	static int times2(int x){
		return x*2;
	}
	static void Main(){
		System.Console.Write($"x*2 = {times2(x)}\n");
	}
}
