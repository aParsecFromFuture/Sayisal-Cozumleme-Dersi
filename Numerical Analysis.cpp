#include <stdio.h>
#include <math.h>

#define ABS(x) ((x < 0)? -(x) : (x))
#define PI 3.14159265358979323846
#define H 0.0001f

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

void fixed_point(float (*g)(float), float p0, float tol, int max_iter)
{
    int iter = 0;
    float p;

    printf("iter\tp(n)\n\n");
    printf("%d.\t%.8f\n", iter, p0);

    while (iter++ < max_iter) {
        p = g(p0);

        printf("%d.\t%.8f\n", iter, p);

        if (ABS(p - p0) < tol) {
            printf("[COMPLETED]\n");
            return;
        }

        p0 = p;
    }
    printf("[ERROR] Max iteration limit exceeded\n");
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
    float p, fa, fp;

    printf("iter\ta\tb\tp\tf(a)\tf(p)\n\n");

    fa = f(a);

    while (iter++ < max_iter) {
        p = (a + b) / 2;
        fp = f(p);

        printf("%d.\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\n", iter, a, b, p, fa, fp);

        if (fp == 0 || ((b - a) / 2) < tol) {
            printf("[COMPLETED]\n");
            return;
        }

        if (fa * fp > 0) {
            a = p;
            fa = fp;
        }
        else {
            b = p;
        }
    }
    printf("[ERROR] Max iteration limit exceeded\n");
}

/*
 *  M) NEWTON'S METHOD
 *
 * Bir fonksiyonun kokunu bulmak icin kullanilir.
 * Bir sayinin degisik derecelerden koklerini hesaplarken de newton ralphson yontemini kullanabiliriz.
 * 
 * Iterasyon ifademiz asagidaki gibidir.
 * 
 * g(x) = x - (f(x) / f'(x))
 * 
 * Karekokunu bulmak istedigimiz sayiya A dersek
 * 
 * sqrt(A) = x
 * A = x^2
 * x^2 - A = 0 elde ederiz
 * 
 * Ifadede yerine yazarsak
 * 
 * g(x) = x - (x^2 - A) / 2x elde ederiz
 * 
 * r'nci dereceden kok icin ise ayni islemleri tekrarlarsak
 * 
 * A^(1/r) = x
 * A = x^r
 * x^r - A = 0 
 *
 * g(x) = ((r - 1)x^r + A) / (rx^(r - 1)) elde ederiz
 * 
 *
 * UYGULANİSİ
 *
 * p0 noktasi secilir.
 * p = p0 - f(p0) / f'(p0) ifadesi ile p noktasi hesaplanir.
 * Eger |p - p0| tolerans degerinden kucukse algoritma sonlandirilir.
 * Degilse p0 = p yapilarak iterasyonlara devam edilir.
 *
 */

void newton_raphson(float (*f)(float), float p0, float tol, int max_iter)
{
    int iter = 0;
    float p;

    printf("%d.\t%.8f\n", iter, p0);

    while (iter++ < max_iter) {
        p = p0 - f(p0) / ((f(p0 + H) - f(p0)) / H);

        printf("%d.\t%.8f\n", iter, p);

        if (ABS(p - p0) < tol) {
            printf("[COMPLETED]\n");
            return;
        }

        p0 = p;
    }
    printf("[ERROR] Max iteration limit exceeded\n");
}

/*
 *  M) SECANT METHOD
 *
 * Bir fonksiyon üzerinde çizeceğimiz sekant doğrusundan yararlanılarak geliştirilmiş bir yöntemdir.
 * Eğri üzerinden iki nokta seçilerek çizilen doğruya sekant doğrusu denir. (x0, f(x0)) ve (x1, f(x1)).
 * x0 ve x1 noktalarından geçen doğrunun x eksenini kestiği noktaya x2 deriz.
 * x2 ile x1 arasındaki mesafe x1 ile x0 arasındaki mesafeden azdır. Bir yaklaşım söz konusudur.
 * Hesaplama genellestirildiginde soyle bir ifade elde edilir:
 * 
 * p(k + 1) = g(p(k - 1), p(k)) = p(k) - f(p(k)) (p(k) - p(k - 1)) / (f(p(k)) - f(p(k - 1))) for k = 0, 1, ...
 * 
 * Sadelestirirsek
 * 
 * p(k + 1) = (p(k - 1) * f(p(k)) - p(k) * f(p(k - 1))) / f(p(k)) - f(p(k - 1))
 * 
 * 
 * UYGULANİSİ
 * 
 * p0 ve p1 noktalari secilir.
 * p = p1 - f(p1) * (p1 - p0) / (f(p1) - f(p0)) ifadesinden p noktasi bulunur.
 * Eğer |p - p1| tolerans değerinden kucukse algoritma sonlandirilir.
 * Degilse p0 = p1, p1 = p yapilarak iterasyonlara devam edilir.
 */

void secant_method(float(*f)(float), float p0, float p1, float tol, int max_iter) {
    int iter = 0;
    float p, q0, q1;

    printf("iter\tp(n)\n\n");

    q0 = f(p0);
    q1 = f(p1);

    while (iter++ < max_iter) {
        p = p1 - q1 * (p1 - p0) / (q1 - q0);

        printf("%d.\t%.8f\n", iter, p);

        if (ABS(p - p1) < tol) {
            printf("[COMPLETED]\n");
            return;
        }

        p0 = p1;
        q0 = q1;
        p1 = p;
        q1 = f(p);
    }
    printf("[ERROR] Max iteration limit exceeded\n");
}

/*
 *  M) HALLEY'S METHOD
 *
 * Newton raphson yonteminin genisletimi ve hizlandirimi olarak alinabilir. 
 * Kokun kat sayisi arttikca koke ulasim hizi azalir.
 * Halley metodu g(x) = x - f(x) / f'(x) ifadesindeki f(x) / f'(x) kısmına bir katsayisi ekler.
 * Bu katsayi birden buyuk oldugu icin koke daha fazla yaklasim elde edilmis olur.
 * 
 * p(k + 1) = p(k) - (f(p(k)) / f'(p(k))) * (1 - (f(p(k)) * f''(p(k))) / (2 * f'(p(k))^2))^(-1) for k = 0, 1, ...
 *
 *
 * UYGULANİSİ
 *
 * p0 noktasi secilir.
 * p = p0 - (f(p0) / f'(p0)) * (1 - (f(p0) * f''(p0)) / (2 * f'(p0)^2))^(-1) ifadesi ile p noktasi hesaplanir.
 * Eger |p - p0| tolerans degerinden kucukse algoritma sonlandirilir.
 * Degilse p0 = p yapilarak iterasyonlara devam edilir.
 * 
 */

void halley_method(float (*f)(float), float p0, float tol, int max_iter)
{
    int iter = 0;
    float p, fp, dfp, ddfp;

    printf("%d.\t%.8f\n", iter, p0);

    while (iter++ < max_iter) {
        fp = f(p0);
        dfp = (f(p0 + H) - f(p0)) / H;
        ddfp = (2 * f(p0 + H) - f(p0) - f(p0 + 2 * H)) / H;
        p = p0 - (fp / dfp) / (1 - (fp * ddfp) / (2*dfp*dfp));

        printf("%d.\t%.8f\n", iter, p);

        if (ABS(p - p0) < tol) {
            printf("[COMPLETED]\n");
            return;
        }

        p0 = p;
    }
    printf("[ERROR] Max iteration limit exceeded\n");
}

float f_1(float x)
{
    return (x * x * x) + 4 * (x * x) - 10;
}

float f_2(float x)
{
    return pow(3 + x - 2 * (x * x), 0.25);
}

float f_3(float x)
{
    return cos(x) - x;
}

float f_4(float x)
{
    return -4 + 17 * (x)-16 * (x * x) + 4 * (x * x * x);
}

int main(int argc, char** argv)
{
    float a, b, tol, max_iter, p0, p1;
    printf("\nf(x) = x^3 + 4x^2 - 10 fonksiyonu icin [-2, 1.5] araliginda bisection metodunun uygulanisi\n\n");
    a = -2.0f;
    b = 1.5f;
    tol = 1e-2;
    max_iter = 100;
    bisection(f_1, a, b, tol, max_iter);

    printf("\nf(x) = (3 + x - 2x^2)^(0.25) fonksiyonu icin p0 = 1 baslangic noktasi uzerinden fixed point uygulanisi\n\n");
    p0 = 1;
    tol = 1e-2;
    max_iter = 100;
    fixed_point(f_2, p0, tol, max_iter);

    printf("\nf(x) = cos(x) - x fonksiyonu icin p0 = pi/4 baslangic noktasi icin newton raphson yontemi\n");
    p0 = PI / 4.0f;
    tol = 1e-6;
    max_iter = 100;
    newton_raphson(f_3, p0, tol, max_iter);

    printf("\nf(x) = cos(x) - x fonksiyonu icin p0 = pi/4 baslangic noktasi icin halley yontemi\n");
    p0 = PI / 4.0f;
    tol = 1e-6;
    max_iter = 100;
    halley_method(f_3, p0, tol, max_iter);

    printf("\nf(x) = -4 + 17x - 16x^2 + 4x^3 fonksiyonu icin p0 = 3, p1 = 2.8 baslangic noktalari icin secant yontemi\n\n");
    p0 = 3.0f;
    p1 = 2.8f;
    tol = 1e-4;
    max_iter = 100;
    secant_method(f_4, p0, p1, tol, max_iter);

    return 0;
}