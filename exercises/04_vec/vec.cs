using System;

public class vec{
	public double x, y, z;
	//constructors:
	public vec(){ x=y=z=0; }
	public vec(double x,double y,double z){ this.x=x; this.y=y; this.z=z;}
	//operators:
	public static vec operator*(vec v, double c){return new vec(c*v.x,c*v.y,c*v.z);}
	public static vec operator*(double c, vec v){return v*c;}
	public static vec operator+(vec u, vec v){
	return new vec(u.x+v.x, u.y+v.y, u.z+v.z);}
	public static vec operator-(vec u){return new vec(-u.x,-u.y,-u.z);}
	public static vec operator-(vec u, vec v){
	return new vec(u.x-v.x, u.y-v.y, u.z-v.z);}
	//methods:
	public void print(string s){Console.Write(s);Console.WriteLine($"{x} {y} {z}");}
	public void print(){this.print("");}
	public double dot(vec other) // to be called as u.dot(v)
	{return this.x*other.x+this.y*other.y+this.z*other.z;}
	public static double dot(vec v,vec w) // to be called as vec.dot(u,v)
	{return v.x*w.x+v.y*w.y+v.z*w.z;}
	public vec cross(vec other) // to be called as u.cross(v)
	{return new vec(this.y*other.z-this.z*other.y,
			this.z*other.x-this.x*other.z,
			this.x*other.y-this.y*other.x);
	}
	public static vec cross(vec u, vec v) // to be called as vec.cross(u, v)
	{return new vec(u.y*v.z-u.z*v.y,
			u.z*v.x-u.x*v.z,
			u.x*v.y-u.y*v.x);
	}
	public double norm(vec other) // to be called as u.norm(v)
	{return Math.Sqrt(Math.Pow(this.x-other.x,2)
		       	+ Math.Pow(this.y-other.y,2)
		       	+ Math.Pow(this.z-other.z,2)); 
	}
	public static double norm(vec u, vec v) // to be called as vec.norm(u, v)
	{return Math.Sqrt(Math.Pow(u.x-v.x,2)
		       	+ Math.Pow(u.y-v.y,2)
		       	+ Math.Pow(u.z-v.z,2)); 
	}
	public static bool approx(double x, double y,double acc = 1e-9, double eps = 1e-9)
	{
		if(Math.Abs(x-y) > acc) return false;
		if(Math.Abs(x-y)/(Math.Abs(x) + Math.Abs(y)) > eps) return false;
		return true;
	}
	public static bool approx(vec u, vec v, double acc = 1e-9, double eps = 1e-9)
	{
		if(!approx(u.x, v.x, acc, eps)) return false;
		if(!approx(u.y, v.y, acc, eps)) return false;
		if(!approx(u.z, v.z, acc, eps)) return false;
		return true;
	}
	public bool approx(vec other, double acc = 1e-9, double eps = 1e-9)
	{return approx(this, other, acc, eps);}
	public override string ToString()
	{return $"{x} {y} {z}";}
}
