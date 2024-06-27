using System;

class main{
static vec[] RandomVectorArray(int N=10){
	Random rand = new Random();
	vec[] VecArray = new vec[N];
	for(int i = 0; i<N; i++){
		double x = (rand.NextDouble() - 0.5) * 10;
		double y = (rand.NextDouble() - 0.5) * 10;
		double z = (rand.NextDouble() - 0.5) * 10;
		VecArray[i] = new vec(x, y, z);
	}
	return VecArray;
}
static void Main(){
	bool passed;
	int FailedTests = 0;
	int N = 10;
	vec ZeroVector = new vec(0.0, 0.0, 0.0);
	vec[] Avecs = RandomVectorArray(N);
	vec[] Bvecs = RandomVectorArray(N);
	/* for(int i = 0; i<N; i++){
	Console.WriteLine(Avecs[i]);
	} */

	Console.WriteLine("Testing A-B = -(B-A) ...");
	passed = true;
	for(int i = 0; i<N; i++){
		if(!vec.approx(Avecs[i] - Bvecs[i], -(Bvecs[i] - Avecs[i]))){
			passed = false;
			FailedTests++;
			break;
		}
	}
	Console.WriteLine(passed ? " ... passed" : " ... failed");

	Console.WriteLine("Testing 2*A = A+A ...");
	passed = true;
	for(int i = 0; i<N; i++){
		if(!vec.approx(2 * Avecs[i], Avecs[i] + Avecs[i])){
			passed = false;
			FailedTests++;
			break;
		}
	}
	Console.WriteLine(passed ? " ... passed" : " ... failed");

	Console.WriteLine("|A|² = A·A ...");
	passed = true;
	for(int i = 0; i<N; i++){
		if(!vec.approx(Math.Pow(Avecs[i].norm(ZeroVector), 2), Avecs[i].dot(Avecs[i]))){
			passed = false;
			FailedTests++;
			break;
		}
	}
	Console.WriteLine(passed ? " ... passed" : " ... failed");

	Console.WriteLine("Testing A⨯B = -B⨯A ...");
	passed = true;
	for(int i = 0; i<N; i++){
		if(!vec.approx(Avecs[i].cross(Bvecs[i]), - (Bvecs[i].cross(Avecs[i])))){
			passed = false;
			FailedTests++;
			break;
		}
	}
	Console.WriteLine(passed ? " ... passed" : " ... failed");

	Console.WriteLine("Testing |A⨯B|² = (A·A)(B·B) - (A·B)² ...");
	passed = true;
	for(int i = 0; i<N; i++){
		if(!vec.approx(
					Math.Pow(vec.norm(Avecs[i].cross(Bvecs[i]), ZeroVector), 2),
					vec.dot(Avecs[i], Avecs[i])*vec.dot(Bvecs[i], Bvecs[i])
					- Math.Pow(Avecs[i].dot(Bvecs[i]), 2))){
			passed = false;
			FailedTests++;
			break;
		}
	}
	Console.WriteLine(passed ? " ... passed" : " ... failed");

	if (FailedTests==0){Console.WriteLine("\nPassed all tests.");}
	else {Console.WriteLine($"\n Failed {FailedTests} tests");}
}
}
