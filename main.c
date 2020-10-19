#include <stdio.h>
#include <math.h>

#define ABS(x) ((x < 0)? -(x) : (x))

/*                      SAYISAL COZUMLEME DERS NOTLARİ
 * 
 *  A) SAYİSAL COZUMLEME METODLARİNİN SİNİFLANDİRİLMASİ
 *  
 *  1.Dogrusal cebir konulari
 *      1.1. Dogrusal denklem sistemlerinin cozumleri
 *      1.2. Ozdeger/ozvektor hesabi
 *      1.3. Matrislerin carpimlarina ayrilmasi
 *  2. Calculus konulari
 *      2.1. Sayisal turev/integral hesabi
 *      2.2. Interpolasyon
 *      2.3. Dogrusal olmayan denklem sistemlerinin cozumleri
 *  3. Istatiki hesaplama konulari
 *      3.1. Polinom yaklasimi
 *      3.2. Egri uydurma
 * 
 *  B) SAYİSAL VE ANALİTİK YONTEMLERİN KARSİLASTİRİLMASİ
 *
 * Muhendislikte cozume ihtiyac duydugumuz problemlerin cogu sayisal yontemlerle ancak cozulebilir.
 * Sayisal yontemleri analitik yontemlerin uygulandigi problemlere uygulamakta mümkün.
 * Analitik yontemler bize kesin degerler uretiyor olsa da sayisal yontemlerle elde edecegimiz cozumler, 
 * ihtiyac duyacagimiz dogruluk icin yeterli olabiliyor. 
 * 
 *  C) DİRECT VE ITERATIVE METOD ARASİNDAKİ FARKLAR
 * 
 *  Direct method
 * 
 * Belirli adimlar icerisinden dogrudan hesaplama yapar. 
 * Ozellikle dogrusal denklem sistemlerinin cozumunu belirli adimlarda elde edersiniz.
 * 
 *  Iterative method
 *  
 * Uygun kosullar altinda dogru cozume ulastiran yontemlerdir. Cozume belirli bir yaklasimla ulasabiliyoruz.
 * Burada iki konu gundeme geliyor. Metodun bizi cozume yaklastirip yaklastirmayacagi ve hangi asamada durdurulacagi.
 * Birincisinin analizi bazi metodlar icin yapilabiliyor, digerleri icinse islemler
 * yapilirken hesaplanan degerlerdeki degisimlere bakarak bir yakinsak cozume ulasilabilirligi belirlenir.
 * 
 *  D) SONLANDİRMA KOSULLARİ
 * 
 * Sonlandirma kosullari icin uc yol kullaniyoruz. 
 * Ornegin x* bir fonksiyon yada polinomun bir koku oldugunu dusunun.
 * xk ise iterative yontemlerden biriyle k'nci iterasyonda elde edilen bir deger oldugunu dusunun. x* fonksiyonun koku olduguundan
 * 0 degerini vermesi beklenir. Ancak xk 0 degerini degil 0'a yakin bir deger verecektir. Bu 0'a yakinlik istedigimiz bir dogrulaga
 * sahipse xk'yi polinomun yaklasik koku olarak alabiliyoruz. k'nci adimda istenilen dogrulaga yada toleransa ulasilip ulasilmadigi
 * degerlendiriliyor. Istenilen dogruluga ulasildiysa algoritma sonlandirilabiliyor. 
 * 
 *  E) HATA TURLERİ
 * 
 *  1.Kesme hatasi
 * 
 * Sonsuz bir serinin toplami yaklasik olarak bu terimlerin sonlu sayisi kullanilarak hesaplandiginda yaklasik bir deger elde edilmis
 * olacak ve burada kesme hatasi meydana gelir. 
 * 
 *  2.Yuvarlama hatasi
 * 
 * Ozellikle bilgisayar sistemlerinde reel sayilarin temsilinde ortaya cikar. Reel sayilar yaklasik
 * olarak bilgisayar ortamlarinda temsil edilirler. Float turu noktadan snora 7, double turu 15 haneye kadar
 * dogru olarak temsil eder.
 *
 *  F) HATA OLCUMU SİNİFLARİ
 * 
 *  Mutlak hata
 * 
 * Iterasyonlar degeri ile gercek deger arasindaki uzakliga mutlak hata denir. Iki ardisik iterasyon degerlerinin farkinin
 * mutlak degerleri alinarak da belirleniyor.
 * 
 *  Bagil hata
 * 
 * Mutlak hataya sonuncu iterasyonda hesaplanan p degeri payda verisi olarak eklenerek bagil hata hesaplaniyor.
 * 
 *  G) TUREVIN LIMIT UZERINDEN GOSTERIMI
 * 
 *  lim      f(x) - f(x0)                lim      f(x0 + h) - f(x0)
 * x -> x0   ------------      veya     h -> 0    ----------------
 *              x - x0                                    h
 *  H) TUREVLENEBILIR FONKSIYONLAR
 * 
 *  1.Ortalama Deger Teoremi
 * 
 * a ve b noktalari arasinda turevlenebilir bir fonksiyon dusunelim. a ve b arasinda oyle bir nokta (c) vardir ki,
 * o noktanin turevi (tanjant dogrusu) a ve b arasinda cizilen sekant dogrusuna paraleldir.
 * 
 *  I) POLINOMLARIN DEGERLENDIRILMESI
 * 
 * p(x) = a0 + a1^x + a2x^2 + ... + anx^n
 * 
 * 1.Yontem (n toplama, (n(n+1) / 2) carpma islemi)
 * 
 * poly = a0
 * for j = 1:n
 *  poly = poly + ajx^j
 * end
 * 
 * 2.Yontem (n toplama, (2n - 1) carpma islemi)
 * 
 * poly = a0 + a1x
 * power = x
 * for j = 2:n
 *  power = x * power
 *  poly = poly + aj * power
 * end
 * 
 * 3.Yontem (n toplama, n carpma islemi)
 * 
 * poly = an
 * for j = n - 1:-1:0
 *  poly = aj + x * poly
 * end
 * 
 */

/*
 *  K) SABIT NOKTA ITERASYONU
 * Sabit noktalari kullanarak gelistirilmis kok bulma yontemidir. Sabit noktalar ile kok hesabina basit nokta iterasyonu da denir.
 * Biz burada sabit nokta iterasyonunu kok bulma iterasyonuna donusturmeye ihtiyac duyariz.
 * Sabit noktanin hesabinda iterasyonlar bizi sabit noktaya tasir. 
 * 
 * UYGULANİSİ
 * 
 * Baslangic olarak p0 noktasi secilir. 
 * p1 = f(p0) ile p1 degeri bulunur.
 * p2 = f(p1) ile p2 degeri bulunur.
 * ...
 * p(k+1) = p(k) ile p(k+1) degeri bulunur.
 * 
 * p(n) - p(n - 1) arasindaki fark tolerans degerinden kucuk olana degin iterasyonlara devam edilir.
 * 
 * NOT! Sabit nokta iterasyonunu yapabilmemiz icin iterasyon araliginda fonksiyoonunun turevi mutlak degerce birin altinda olmali.
 * NOT! Bir fonksiyonun y = x dogrusu ile kesisim noktalarina o fonksiyonun sabit noktalari diyoruz.
 */

void fixed_point(float (*g)(float), float x, float tol, int max_iter)
{
    int iter = 0;
    float err, new_x;
    
    printf("p0:\t%.5f\n", x);
    
    do{
        new_x = g(x);
        err = ABS(x - new_x);
        x = new_x;
        printf("p%d:\t%.5f\n", iter + 1, x);
    }while(iter++ < max_iter && err >= tol);
    
    printf("[COMPLETED]\n");
}

/*
 *  L) BISECTION YONTEMI
 * 
 * Eger bir fonksiyonu verilen bir a - b araliginda bir koku varsa, fonksiyonun a noktasinda aldigi deger ile b
 * noktasinda aldigi degerin zit isaretli olmasi gerekir. Buna "orta deger teoremi" denir.
 * Kokun varligini anlayabilmemiz icin iterasyonlara baslarken bu kosula bakariz. 
 * Yalniz bu kosulun saglamamasi kok olmadigi anlamina gelmez. 
 * Birden fazla kok oldugu anlamina da gelebilir. Ancak en az bir tane kok olduguna bu kosula bakarak emin olabiliriz. 
 * 
 * UYGULANİSİ
 * 
 * a - b araliginin orta noktasi c'yi hesaplariz.
 * Eger f(a) * f(c) sifirdan kucukse b = c olur. 
 * Eger f(b) * f(c) sifirdan kucukse a = c olur. 
 * Ikı kosuldan biri saglanmaz ise (f(c) = 0 olursa) algoritma sonlandirilir. 
 * a - b araliginin yarisi tolerans degerinden kucuk olana kadar iterasyonlara devam edilir.
 * 
 */

void bisection(float (*f)(float), float a, float b, float tol, int max_iter)
{
    int iter = 0;
    float c, err, f_a, f_b, f_c;
    
    f_a = f(a);
    f_b = f(b);
    
    if(f_a * f_b >= 0){
        printf("[ERROR] f(a) * f(b) should be negative!\n");
        return;
    }
    
    err = ABS(a - b) * 0.5f;
    
    printf("\ta\t\tb\t\tc\t\tf(a)\t\tf(b)\t\tf(c)\n\n");
    
    while(iter++ < max_iter && err >= tol){
        c = (a + b) * 0.5f;
        f_c = f(c);
        
        printf("%d.\t%.4f\t\t%.4f\t\t%.4f\t\t%.4f\t\t%.4f\t\t%.4f\n", iter, a, b, c, f_a, f_b, f_c);
        
        if(f_a * f_c < 0){
            b = c;
            f_b = f_c;
        }else if(f_b * f_c < 0){
            a = c;
            f_a = f_c;
        }else{
            printf("[ROOT FOUND] %.2f", c);
            return;
        }
        err = ABS(a - b) * 0.5f;
    }
    printf("[COMPLETED]\n");
}

float f_1(float x)
{
    return (x * x * x) + 4 * (x * x) - 10;
}

float f_2(float x)
{
    return pow(3 + x - 2 * (x * x), 0.25);
}

int main(int argc, char **argv)
{
    printf("\nf(x) = x^3 + 4x^2 - 10 fonksiyonu icin [-2, 1.5] araliginda 10^-2 hata payi ile bisection metodunun uygulanisi\n\n");
    float a = -2.0f;
    float b = 1.5f;
    float tol = 1e-2;
    int max_iter = 100;
	bisection(f_1, a, b, tol, max_iter);
    
    
    printf("\nf(x) = (3 + x - 2x^2)^(0.25) fonksiyonu icin p0 = 1 baslangic noktasi uzerinden 10^-2 hata payi ile fixed point uygulanisi\n\n");
    float p0 = 1;
    fixed_point(f_2, p0, tol, max_iter);
    
	return 0;
}
