using System;
using System.Collections;
using System.Collections.Generic;

public class GenericList<T> : IEnumerable<T>{
	private T[] items;
	private int size;

	public GenericList(){items = new T[0]; size = 0;}
	public void Add(T item){
		T[] newList = new T[size+1];
		System.Array.Copy(items, newList, size);
		newList[size] = item;
		items = newList;
		size++;
	}
	public IEnumerator<T> GetEnumerator(){
		for (int i = 0; i<size; i++){
			yield return items[i];
		}
	}

	IEnumerator IEnumerable.GetEnumerator() {
		return GetEnumerator();
	}
}

public class Program{
public static void Main(){
	var list = new GenericList<double>();
	char[] delimiters = {' ', '\t'};
	var split_options = StringSplitOptions.RemoveEmptyEntries;

	for (string line = Console.ReadLine(); line!= null; line = Console.ReadLine()){
		var numbers = line.Split(delimiters, split_options);

		foreach(var number in numbers){
			list.Add(double.Parse(number));
		}
	}

	foreach (var item in list){
		Console.WriteLine($"Here is a number: {item:e4}");
	}
} // Main
} // Program
