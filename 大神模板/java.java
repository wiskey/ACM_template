方法 
valueOf 
 public static BigInteger valueOf(long val)返回一个是指定值的 BigInteger 。 给该工厂提供首选的 (long) 构造子是因为它允许对频繁使用的 BigIntegers (如 0 和 1)的重用，这种操作不需要输出常数。 

add 
 public BigInteger add(BigInteger val) throws ArithmeticException
返回一个 BigInteger ，其值是 (this + val) 。 
subtract 

 public BigInteger subtract(BigInteger val)返回一个 BigInteger ，其值是 (this - val) 。 
multiply 

 public BigInteger multiply(BigInteger val)返回一个 BigInteger ，其值是 (this * val) 。 
divide 

 public BigInteger divide(BigInteger val) throws ArithmeticException
返回一个 BigInteger ，其值是 (this / val) 。 如果 val == 0 ，则抛出 ArithmeticException 。 
remainder 

 public BigInteger remainder(BigInteger val) throws ArithmeticException
返回一个 BigInteger ，其值是 (this % val) 。 如果 val == 0，则抛出 ArithmeticException 。 
divideAndRemainder 

 public BigInteger[] divideAndRemainder(BigInteger val) throws ArithmeticException
返回一个包含两个 BigIntegers 的数组。 返回值的第一个 ([0]) 元素是商 (this / val), 第二个 ([1]) 元素是余数 (this % val) 。如果 val == 0 ，则抛出 ArithmeticException 。 
pow 

 public BigInteger pow(int exponent) throws ArithmeticException
返回一个 BigInteger ，其值是 (this ** exponent) 。 如果 exponent <0(因为该操作将产生一个非整型值)，则抛出 arithmeticexception 。注意：指数是一个整型数而不是 biginteger 。 
gcd 

 public BigInteger gcd(BigInteger val)返回值为 abs(this) 和 abs(val) 最大公分母的 BigInteger 。 如果 this == 0 && val == 0，则返回 0。 
abs 

 public BigInteger abs()返回一个 BigInteger ，它是该数值的绝对值。 
negate 

 public BigInteger negate()返回一个 BigInteger ，其值是 (-1 * this ) 。 
signum 

 public int signum()返回该数值的符号数 (根据该数的值是正、零或负返回 -1 、 0 或 1 ) 。 
mod 

 public BigInteger mod(BigInteger m)返回一个 BigInteger ，其值是 this mod m 。 如果 m <= 0，则抛出 arithmeticexception。 
modPow 

 public BigInteger modPow(BigInteger exponent,                          BigInteger m)返回一个 BigInteger ，其值是 (this ** exponent) mod m 。 如果 exponent == 1, 返回值是 (this mod m) 。 如果 exponent <0 ， 返回值是 (this ** -exponent)的模多重逆。如果 m <="0，则抛出" arithmeticexception 。 

modInverse 

 public BigInteger modInverse(BigInteger m) throws ArithmeticException
返回 this 取模 m 的模多重逆。 如果 m<= 0 或 this 没有多重逆 mod m (比如 gcd(this, m) !="1)，则抛出" arithmeticexception 。 
shiftLeft 

 public BigInteger shiftLeft(int n)返回一个 BigInteger ，其值是 (this << n)。计算 floor(this * 2**n)。 
shiftRight 

 public BigInteger shiftRight(int n)返回一个 BigInteger ，其值是 (this >> n)。 执行符号扩展(计算 floor(this / 2**n))。 
and 

 public BigInteger and(BigInteger val)返回一个 BigInteger ，其值是 (this & val) 。 (如果 this 和 val 二者都是负的，则该方法返回一个负数。) 
or 

 public BigInteger or(BigInteger val)返回一个 BigInteger ，其值是 (this | val) 。 (如果 this 和 val 二者之一是负的，则该方法返回一个负数。) 
xor 

 public BigInteger xor(BigInteger val)返回一个 BigInteger ，其值是 (this ^ val) 。 (如果 this 和 val 二者只有一个是负的，则该方法返回一个负数。) 
not 

 public BigInteger not()返回一个 BigInteger ，其值是 (~this) 。 (如果 this 是非负的，则该方法返回一个负数。) 
andNot 

 public BigInteger andNot(BigInteger val)返回一个 BigInteger ，其值是 (this & ~val) 。 该方法等价于 and(val.not())，它是进行掩模操作的便捷方法。 (如果 this 是负数并且 val 是正数，则该方法返回一个负数。) 
testBit 

 public boolean testBit(int n) throws ArithmeticException
如果设置了指定位则返回 true 。 (计算 ((this & (1 << n)) !=0))。如果 n < 0，则抛出 arithmeticexception。 
setBit 

 public BigInteger setBit(int n) throws ArithmeticException
返回一个 BigInteger ，其值等于该数被设置指定位后所得值(计算 (this | (1 << n)))。 如果 n < 0 ，则抛出 arithmeticexception 。 
clearBit 

 public BigInteger clearBit(int n) throws ArithmeticException
返回一个 BigInteger ，其值等于该数指定位清零后所得值(计算 (this & ~(1 << n)))。 如果 n < 0，则抛出 arithmeticexception 。 
flipBit 

 public BigInteger flipBit(int n) throws ArithmeticException
返回一个 BigInteger ，其值等于该数指定位取反后所得值(计算 (this ^ (1 << n)))。 如果 n < 0，则抛出 arithmeticexception 。 
getLowestSetBit 

 public int getLowestSetBit()返回该数最右端 (lowest-order)是 1 的位的索引 (即就是距最右端 1 位间的 0 位的个数 ) 。 如果该数没有 1 位，则返回 -1 (计算 (this==0? -1 : log2(this & -this)))。 
bitLength 

 public int bitLength()返回该数的最小二进制补码表示的位的个数， 即 *不包括* 符号位 (ceil(log2(this <0 ? -this : this + 1)))。对正数来说，这等价于普通二进制表示的位的个数。 
bitCount 

 public int bitCount()返回该数的二进制补码表示中不包扩符号位在内的位的个数。 该方法在 BigIntegers 之上实现位向量风格的集合时很有用。 
isProbablePrime 

 public boolean isProbablePrime(int certainty)如果该 BigInteger 可能是素数，则返回 true ；如果它很明确是一个合数，则返回 false 。 参数 certainty 是对调用者愿意忍受的不确定性的度量：如果该数是素数的概率超过了 1 - 1/2**certainty方法，则该方法返回 true 。执行时间正比于参数确定性的值。 
compareTo 

 public int compareTo(BigInteger val)根据该数值是小于、等于、或大于 val 返回 -1、0 或 1 。该方法在六个逻辑比较运算符 (<, ="=,">, >=, !=, <= ) 的操作中作为首选方法。进行这些比较的建议方法是：(x.compareto(y) 0) ，其中 是六个比较符中的一个。 
equals 

 public boolean equals(Object x)如果 x 是一个 BigInteger 并且等于该数则返回 true 。 提供该方法的目的是使 BigIntegers 可以作为散列码关键字使用。 
覆盖： 
类 Object 中的 equals 
min 

 public BigInteger min(BigInteger val)返回 BigInteger ，其值是 this 和 val 中的较小者。 若值相同，则两者都可能被返回。 
max 

 public BigInteger max(BigInteger val)返回 BigInteger ，其值是 this 和 val 中的较大者。 若值相同，则两者都可能被返回。 
hashCode 

 public int hashCode()为该对象计算一个散列码。 
覆盖： 
类 Object 中的 hashCode 
toString 

 public String toString(int radix)返回表示该数的字符串，基数为给定基数 。 如果基数超出了 Character.MIN_RADIX(2) 到 Character.MAX_RADIX(36)(包括两者在内) 的范围，它会缺省设定为 10 (类似于 Integer.toString 的情形)。使用由 Character.forDigit 提供的数字到字符的映射，并且如果适当的话，还可以前置一个负号。 该表示法同 (String, int) 构造子兼容。 
toString 

 public String toString()返回表示该数的字符串，基数为 10 。 使用由 Character.forDigit 提供的数字到字符的映射，并且如果合适的话，还可以前置一个负号。 该表示法同 (String) 构造子兼容，并且允许用 Java 的 + 运算符做字符串连接操作。 
覆盖： 
类 Object 中的 toString 
toByteArray 

 public byte[] toByteArray()返回该数值的二进制补码表示。 数组是 big-endian 格式(即最有效字节在位置 [0]) 。该数组包含了需要表示该数的最小的字节数 (ceil((this.bitLength() + 1)/8)) 。该表示法同(byte[])构造子兼容。 
intValue 

 public int intValue()把该数转换为 int 值。 标准的限制原语转换同《Java 语言规范》所定义的一样。 
覆盖： 
类 Number 中的 intValue 
longValue 

 public long longValue()把该数转换为 long 型。标准的限制原语转换同《Java 语言规范》定义的一样。 
覆盖： 
类 Number 中的 longValue 
floatValue 

 public float floatValue()把该数转换为 float 型。该操作类似于《Java 语言规范》中定义的 double-to-float 限制原语转换：如果数值太大以致不能表示为一个浮点数时，它将被转换为合适的无穷大或负无穷大。 
覆盖： 
类 Number 中的 floatValue 
doubleValue 

 public double doubleValue()把该数转换为 double 型。 该操作类似于《Java 语言规范》中定义的 double-to-float 的限制原语转换：如果数值太大以致不能表示为一个双精度数，它将被转换为适当的无穷大或负无穷大。 

import java.math.BigInteger;
public class BigIntegerGet {
    public String getAdd(String Str1,String Str2){
        String Str3=new String();
        BigInteger BigInt1=new BigInteger(Str1);
        BigInteger BigInt2=new BigInteger(Str2);
        BigInt1=BigInt1.add(BigInt2);
        Str3=BigInt1.toString();
        return Str3;
    }
    public String getSubtract(String Str1,String Str2){
        String Str3=new String();
        BigInteger BigInt1=new BigInteger(Str1);
        BigInteger BigInt2=new BigInteger(Str2);
        BigInt1=BigInt1.subtract(BigInt2);
        Str3=BigInt1.toString();       
        return Str3;
    }
    public String getMultiply(String Str1,String Str2){
        String Str3=new String();
        BigInteger BigInt1=new BigInteger(Str1);
        BigInteger BigInt2=new BigInteger(Str2);
        BigInt1=BigInt1.multiply(BigInt2);
        Str3=BigInt1.toString();       
        return Str3;
    }
    public String getDivide(String Str1,String Str2){
        String Str3=new String();
        BigInteger BigInt1=new BigInteger(Str1);
        BigInteger BigInt2=new BigInteger(Str2);
        BigInt1=BigInt1.divide(BigInt2);
        Str3=BigInt1.toString();       
        return Str3;
    }
    public String getRemainder(String Str1,String Str2){//％
        String Str3=new String();
        BigInteger BigInt1=new BigInteger(Str1);
        BigInteger BigInt2=new BigInteger(Str2);
        BigInt1=BigInt1.remainder(BigInt2);
        Str3=BigInt1.toString();       
        return Str3;
    }
    public String getGcd(String Str1,String Str2){
        String Str3=new String();
        BigInteger BigInt1=new BigInteger(Str1);
        BigInteger BigInt2=new BigInteger(Str2);
        BigInt1=BigInt1.gcd(BigInt2);
        Str3=BigInt1.toString();       
        return Str3;
    }
    public String getPow(String Str1,String Str2){
        String Str3=new String();
        BigInteger BigInt1=new BigInteger(Str1);
        int Int2=Integer.valueOf(Str2);
        BigInt1=BigInt1.pow(Int2);
        Str3=BigInt1.toString();       
        return Str3;
    }
    public String getMod(String Str1,String Str2){
        String Str3=new String();
        BigInteger BigInt1=new BigInteger(Str1);
        BigInteger BigInt2=new BigInteger(Str2);
        BigInt1=BigInt1.mod(BigInt2);
        Str3=BigInt1.toString();        
        return Str3;
    }
    public String getModInverse(String Str1,String Str2){
        String Str3=new String();
        BigInteger BigInt1=new BigInteger(Str1);
        BigInteger BigInt2=new BigInteger(Str2);
        BigInt1=BigInt1.modInverse(BigInt2);
        Str3=BigInt1.toString();        
        return Str3;
    }
    public String getMax(String Str1,String Str2){
        String Str3=new String();
        BigInteger BigInt1=new BigInteger(Str1);
        BigInteger BigInt2=new BigInteger(Str2);
        BigInt1=BigInt1.max(BigInt2);
        Str3=BigInt1.toString();        
        return Str3;
    }
    public String getMin(String Str1,String Str2){
        String Str3=new String();
        BigInteger BigInt1=new BigInteger(Str1);
        BigInteger BigInt2=new BigInteger(Str2);
        BigInt1=BigInt1.min(BigInt2);
        Str3=BigInt1.toString();        
        return Str3;
    }
    public int getHashcode(String Str){
        int hash=-1;
        BigInteger BigInt=new BigInteger(Str);
        hash=BigInt.hashCode();
        return hash;
    }
    public boolean getIsProbablePrime(String Str,int certainty){//是素数概率为1 - 1/2certainty
        boolean flag=false;
        BigInteger BigInt=new BigInteger(Str);
        flag=BigInt.isProbablePrime(certainty);
        return flag;
    }
}
import java.math.BigDecimal; 

public class Arith { 

/** 
* 由于Java的简单类型不能够精确的对浮点数进行运算，这个工具类提供精 
* 确的浮点数运算，包括加减乘除和四舍五入。 
*/ 
//默认除法运算精度 
private static final int DEF_DIV_SCALE = 10; 
    
//这个类不能实例化 
private Arith(){ 
} 

    /** 
     * 提供精确的加法运算。 
     * @param v1 被加数 
     * @param v2 加数 
     * @return 两个参数的和 
     */ 
    public static double add(double v1,double v2){ 
        BigDecimal b1 = new BigDecimal(Double.toString(v1)); 
        BigDecimal b2 = new BigDecimal(Double.toString(v2)); 
        return b1.add(b2).doubleValue(); 
    } 
    /** 
     * 提供精确的减法运算。 
     * @param v1 被减数 
     * @param v2 减数 
     * @return 两个参数的差 
     */ 
    public static double sub(double v1,double v2){ 
        BigDecimal b1 = new BigDecimal(Double.toString(v1)); 
        BigDecimal b2 = new BigDecimal(Double.toString(v2)); 
        return b1.subtract(b2).doubleValue(); 
    } 
    /** 
     * 提供精确的乘法运算。 
     * @param v1 被乘数 
     * @param v2 乘数 
     * @return 两个参数的积 
     */ 
    public static double mul(double v1,double v2){ 
        BigDecimal b1 = new BigDecimal(Double.toString(v1)); 
        BigDecimal b2 = new BigDecimal(Double.toString(v2)); 
        return b1.multiply(b2).doubleValue(); 
    } 

    /** 
     * 提供（相对）精确的除法运算，当发生除不尽的情况时，精确到 
     * 小数点以后10位，以后的数字四舍五入。 
     * @param v1 被除数 
     * @param v2 除数 
     * @return 两个参数的商 
     */ 
    public static double div(double v1,double v2){ 
        return div(v1,v2,DEF_DIV_SCALE); 
    } 

    /** 
     * 提供（相对）精确的除法运算。当发生除不尽的情况时，由scale参数指 
     * 定精度，以后的数字四舍五入。 
     * @param v1 被除数 
     * @param v2 除数 
     * @param scale 表示表示需要精确到小数点以后几位。 
     * @return 两个参数的商 
     */ 
    public static double div(double v1,double v2,int scale){ 
        if(scale<0){ 
            throw new IllegalArgumentException( 
                "The scale must be a positive integer or zero"); 
        } 
        BigDecimal b1 = new BigDecimal(Double.toString(v1)); 
        BigDecimal b2 = new BigDecimal(Double.toString(v2)); 
        return b1.divide(b2,scale,BigDecimal.ROUND_HALF_UP).doubleValue(); 
    } 

    /** 
     * 提供精确的小数位四舍五入处理。 
     * @param v 需要四舍五入的数字 
     * @param scale 小数点后保留几位 
     * @return 四舍五入后的结果 
     */ 
    public static double round(double v,int scale){ 
        if(scale<0){ 
            throw new IllegalArgumentException( 
                "The scale must be a positive integer or zero"); 
        } 
        BigDecimal b = new BigDecimal(Double.toString(v)); 
        BigDecimal one = new BigDecimal("1"); 
        return b.divide(one,scale,BigDecimal.ROUND_HALF_UP).doubleValue(); 
    } 
    
   /** 
    * 提供精确的类型转换(Float) 
    * @param v 需要被转换的数字 
    * @return 返回转换结果 
    */ 
    public static float convertsToFloat(double v){ 
    BigDecimal b = new BigDecimal(v); 
    return b.floatValue(); 
    } 
    
    /** 
* 提供精确的类型转换(Int)不进行四舍五入 
* @param v 需要被转换的数字 
* @return 返回转换结果 
*/ 
public static int convertsToInt(double v){ 
BigDecimal b = new BigDecimal(v); 
    return b.intValue(); 
} 

/** 
* 提供精确的类型转换(Long) 
* @param v 需要被转换的数字 
* @return 返回转换结果 
*/ 
public static long convertsToLong(double v){ 
BigDecimal b = new BigDecimal(v); 
    return b.longValue(); 
} 

/** 
* 返回两个数中大的一个值 
* @param v1 需要被对比的第一个数 
* @param v2 需要被对比的第二个数 
* @return 返回两个数中大的一个值 
*/ 
public static double returnMax(double v1,double v2){ 
BigDecimal b1 = new BigDecimal(v1); 
BigDecimal b2 = new BigDecimal(v2); 
    return b1.max(b2).doubleValue(); 
} 

/** 
* 返回两个数中小的一个值 
* @param v1 需要被对比的第一个数 
* @param v2 需要被对比的第二个数 
* @return 返回两个数中小的一个值 
*/ 
public static double returnMin(double v1,double v2){ 
BigDecimal b1 = new BigDecimal(v1); 
BigDecimal b2 = new BigDecimal(v2); 
    return b1.min(b2).doubleValue(); 
} 

/** 
* 精确对比两个数字 
* @param v1 需要被对比的第一个数 
* @param v2 需要被对比的第二个数 
* @return 如果两个数一样则返回0，如果第一个数比第二个数大则返回1，反之返回-1 
*/ 
public static int compareTo(double v1,double v2){ 
BigDecimal b1 = new BigDecimal(v1); 
BigDecimal b2 = new BigDecimal(v2); 
    return b1.compareTo(b2); 
} 
} 
高精度乘方
import java.io.*;
import java.util.*;
import java.math.*;

public class Main{
	public static void main(String args[])throws Exception{
		Scanner cin=new Scanner(System.in);
		int n;
		String ans;
		BigDecimal a;
		while(cin.hasNext()){
			a=cin.nextBigDecimal();
			n=cin.nextInt();
			a=a.pow(n).stripTrailingZeros();
			ans=a.toPlainString();
			if(ans.substring(0,2).startsWith("0."))
		    	ans=ans.substring(1);
			System.out.println(ans);
		}
	}
}

1.String (char a[]) 用一个字符数组a 创建一个字符串对象,如 
  char a[]={‘b’,’o’,’y’}; 
  String s=new String(a); 
  上述过程相当于 String s= "boy"; 

2. String(char a[],int startIndex,int count) 提取字符数组a 中的一部分字符创建一个字符串对象 ,参数startIndex 和count 分别指定在a 中提取字符的起始位置和从该位置开始截取的字符个数,例如 
char a[]={‘s’,’t’,’b’,’u’,’s’,’n’}; 
String s=new String(a,2,3); 
上述过程相当于 String s= "bus"; 

3 equals方法 equalsIgnoreCase方法忽略大小写 
字符串对象调用String类中的public boolean equals(String s)方法比较当前字符串对象的实体是否与参数指定的字符串s的实体相同.如 
String tom=new String( "we are students"); 
String boy=new String( "We are students"); 
String jerry= new String("we are students"); 

tom.equals(boy)的值是false, 
tom.equals(jerry)的值是 true. 

4 startsWith,endsWith方法 
String tom= "220302620629021",jerry= "21079670924022"; 
tom.startsWith("220")的值是true 
jerry.startsWith("220")的值是false. 

5 regionMatches方法 
字符串调用 
public boolean regionMatches(int firstStart,String other,int ortherStart,int length) 
方法,从当前字符串参数firstStart指定的位置开始处,取长度为length的一个子串,并将这个子串和参数other 指定的一个子串进行比较,其中,other 指定的子串是从参数othertStart 指定的位置开始,从other中取长度为length的一个子串.如果两个子串相同该方法就返回true,否则返回false. 
使用该方法的重载方法 
public boolean regionMatches(boolean b,int firstStart,String other,int ortherStart,intlength) 
可以通过参数b决定是否忽略大小写,当b取true时,忽略大小写. 


6 compareTo,compareToIgnoreCase方法 
字符串对象可以使用String类中的 
public int compareTo(String s)方法,按辞典序与参数s指定的字符串比较大小. 
如果当前字符串与s相同,该方法返回值0 
如果当前字符串对象大于s,该方法返回正值 
如果小于s,该方法返回负值. 

7 public int indexOf (String s) 字符串调用该方法从当前字符串的始检索字符串s,并返回首次出现s的位置.如果没有检索到字符串s,该方法返回的值是-1. 

public int indexOf(String s ,int startpoint) 字符串调用该方法从当前字符串的startpoint 位置处开始检索字符串s,并返回首次出现s 的位置.如果没有检索到字符串s,该方法返回的值是-1. 

public int lastIndexOf (String s) 字符串调用该方法从当前字符串的头开始检索字符串s,并返回最后出现s的位置.如果没有检索到字符串s,该方法返回的值是-1. 

public int lastIndexOf(String s ,int startpoint) 字符串调用该方法从当前字符串的startpoint 位置处开始检索字符串s,并返回最后出现s 的位置.如果没有检索到字符串s,该方法返回的值是-1. 

8 public String substring(int startpoint) 字符串对象调用该方法获得一个当前字符串的子串,该子串是从当前字符串的startpoint处截取到最后所得到的字符串. 

public String substring(int start ,int end) 字符串对象调用该方法获得一个当前字符串的子串,该子串是从当前字符串的start 处截取到end 处所得到的字符串,但不包括end处所对应的字符. 

9 public String replace(char oldChar,char newChar) 字符串对象s调用该方法可以获得一个串对象,这个串对象是用参数newChar 指定的字符替换s 中由oldChar 指定的所有字符而得到的字符串. 

public String replaceAll(String old ,String new) 字符串对象s调用该方法可以获得一个串对象,这个串对象是通过用参数new指定的字符串替换s 中由old 指定的所有字符串而得到的字符串. 

Public String trim() 一个字符串s 通过调用方法trim()得到一个字符串对象,该字符串对象是s去掉前后空格后的字符串. 

10 转化为整型 
java.lang包中的Integer类调用其类方法 
public static int parseInt(String s) 
可以将 “数字” 格式的字符串,如"12387",转化为int型数据.例如 
int x; 
String s="6542"; 
x=Integer.parseInt("6542"); 
类似地,使用java.lang包中的Byte,Short,Long类调相应的类方法 
public static byte parseByte(String s) 
public static short parseShort(String s) 
public static long parseLong(String s) 
可以将“数字”格式的字符串,转化为相应的基本数据类型 

11 转化为float型或double型java.lang包中的Float类调用其类方法 
  public static int parseFloat (String s)可以将 “数字”格式的字符串,如"12387.8976",转 化为float 型数据.例如 
  float n=Float.parseFloat("12387.8976") 
或 
  String s= new String(“12387.8976”); 
  float n=Float.parseFloat(s) 

12 public StringvalueOf( byte n) 
  public StringvalueOf (int n) 
  public StringvalueOf (long n) 
  public StringvalueOf (float n) 
  public StringvalueOf (double n) 

将形如123,1232.98等数值转化为字符串,如 
String str=String.valueOf(12313.9876); 
float x=123.987f; 
String temp=String.valueOf(x); 

13 将字符串中的字符拷贝到字符数组 
public void getChars(int start,int end,char c[],int offset ) 字符串调用getChars 方法将当前字符串中的一部分字符拷贝到参数c 指定的数组中,将字符串中从位置start 到end-1位置上的字符拷贝的数组c中,并从数组c 的offset处开始存放这些字符.需要注意的是,必须保证数组c能容纳下要被拷贝的字符.

Java中用split函数进行分割字符串。  
1.语法如下 

String.split(sourceStr,maxSplit) 

String.split(sourceStr) 

参数说明：sourceStr是被分割的字符串，maxSplit是最大的分割数 

返回值说明：split函数的返回值是一个字符串数组String[] 

2.示例代码 

package wang48.jiaocheng; 
public class StringSplit  
{ 
 public static void main(String[]args) 
 { 
  String sourceStr="1,2,3,4,5"; 
  String[] sourceStrArray=sourceStr.split(","); 
  for(int i=0;i<sourceStrArray.length;i++) 
  { 
   System.out.println(sourceStrArray[i]); 
  } 
   
  //最多分割出3个字符串 
  int maxSplit=3; 
  sourceStrArray=sourceStr.split(",",maxSplit); 
  for(int i=0;i<sourceStrArray.length;i++) 
  { 
   System.out.println(sourceStrArray[i]); 
  } 
   
 } 

} 

输出结果： 

1 
2 
3 
4 
5 
1 
2 
3,4,5 

StringTokenizer用法
1。StringTokenizer属于java.util类，引用的时候导入类：import java.util.*;2。构造函数：          1）StringTokenizer（String str）：构造一个解析str的对象，这时候默认的Token间隔符有“空格”、“制表符(‘\t’)”、“换行符(‘\n’)”、“回车符(‘\r’)”。        2）StringTokenizer(String str, String delim) ：构造一个用来解析str的StringTokenizer对象，并提供一个指定的分隔符，如“，”等。     3）StringTokenizer(String str, String delim, boolean returnDelims) ：构造一个用来解析str的StringTokenizer对象，并提供一个指定的分隔符（同2），同时，指定是否返回分隔符。3。方法。 
     说明： 所有方法均为public； 
             1）.int countTokens() ：返回nextToken方法被调用的次数。如果采用构造函数1和2，返回的就是分隔符数量。 
        2.）boolean hasMoreTokens() ：返回是否还有分隔符。 
       3.）boolean hasMoreElements() ：结果同上。 
       4.）String nextToken() ：返回从当前位置到下一个分隔符的字符串。 
       5.）Object nextElement() ：结果同4。 
     6.）String nextToken(String delim) ：与4类似，以指定的分隔符返回结果。
下面有两个例子，在网上找得，没有分析，我在这里分析一下。例1： 
代码: 
      String s = new String("The Java platform is the ideal platform for network computing"); 
      StringTokenizer st = new StringTokenizer(s); 
      System.out.println( "Token Total: " + st.countTokens() ); 
      while( st.hasMoreElements() ){ 
         System.out.println( st.nextToken() ); 
　　　　　　　　　　　} 
结果为： 
Token Total: 10 
The 
Java 
platform 
is 
the 
ideal 
platform 
for 
network 
computing 
这个很简单，就用不找大说了，构造函数1和 方法1、和2。例2: 
代码: 
      String s = new String("The=Java=platform=is=the=ideal=platform=for=network=computing"); 
      StringTokenizer st = new StringTokenizer(s,"=",true); //当为false时看看是什么结果？
      System.out.println( "Token Total: " + st.countTokens() ); 
      while( st.hasMoreElements() ){ 
         System.out.println( st.nextToken() ); 
      } 
结果为： 
Token Total: 19 （false时为10）
The 
= （false时就没有了“＝”）
Java 
= 
platform 
= 
is 
= 
the 
= 
ideal 
= 
platform 
= 
for 
= 
network 
= 
Computing

Java的输入结束判断
Scanner in = new Scanner(System.in);
while(in.hasNext()){
System.out.println(in.next());
}

BufferedReader in=new BufferedReader(new InputStreamReader(System.in));
while(true){
     temp=in.readLine();
     if(temp==null) break;
}

Java 文件读入
Scanner in=new Scanner(new File("data.in"));
BufferedReader in = new BufferedReader( new FileReader("data.in") );

Java文件输出
PrintWriter out=new PrintWriter( "test.out" );
out.println(a.add(b));
out.close(); 