# SOL103 vs SOL111 vs SOL112

Tujuan catatan ini adalah membuat alur ketiga solution sequence lebih mudah dipahami, terutama karena `SOL111` dan `SOL112` sama-sama memakai fondasi modal tetapi tujuan akhirnya berbeda dari `SOL103`.

## Inti singkat

- `SOL103` = modal analysis saja
- `SOL111` = response spectrum berbasis modal
- `SOL112` = response spectrum berbasis modal juga, tetapi dipakai di workflow kita untuk deck/dekode yang mengarah ke response spectrum run dan output RS

Dengan kata lain:

- `SOL103` berhenti setelah mode, frekuensi, periode, participation factor, dan effective mass diketahui
- `SOL111/SOL112` memakai hasil modal itu sebagai bahan untuk menghitung respons spektrum

## Gambaran besar

```text
struktur + boundary + mass
        |
        v
   eigen extraction
        |
        v
 mode shapes / freq / period
        |
        +--> SOL103 selesai di sini
        |
        v
 modal participation factor
 effective modal mass
        |
        v
 response per mode dari spectrum
        |
        v
 hasil arah dasar: RSX / RSY / RSZ
        |
        v
 optional combination
   - SRSS
   - CQC
   - orthogonal percentage rule
        |
        v
 final response spectrum output
```

## SOL103

### Makna

`SOL103` adalah modal analysis murni.

### Input yang dipikirkan user

Biasanya user memberi:

- struktur
- massa
- boundary condition
- parameter eigensolver
- jumlah mode

### Keluaran utama

- eigenvalue
- frequency
- period
- mode shape
- modal participation factor
- effective modal mass

### Apa yang tidak dilakukan SOL103

`SOL103` belum menghitung respons spektrum.

Artinya:

- belum ada `RSX`
- belum ada `RSY`
- belum ada `RSZ`
- belum ada `SRSS`
- belum ada `CQC`

### Cara berpikir sederhana

Kalau pertanyaannya adalah:

> struktur ini punya mode apa saja?

maka itu `SOL103`.

## SOL111

### Makna

Dalam konsep yang kita pakai sekarang, `SOL111` adalah **response spectrum analysis berbasis modal**.

Jadi dia bukan direct transient, dan bukan sekadar modal summary. Dia memakai hasil modal untuk menghitung respons spektrum.

### Alur

```text
SOL111
  -> hitung mode / eigenvector
  -> hitung modal participation factor
  -> hitung effective modal mass
  -> baca response spectrum
  -> hitung response tiap mode
  -> susun hasil arah dasar (RSX/RSY/RSZ)
  -> kalau diminta, lakukan combo
  -> tulis hasil final
```

### Kenapa SOL111 tetap punya tabel modal

Karena fondasinya tetap modal.

Itulah sebabnya F06 `SOL111` bisa tetap memuat:

- frequency
- period
- participation factor
- effective modal mass

Tabel itu bukan aneh. Justru itu menjelaskan bahan baku yang dipakai untuk menghitung respons spektrum.

### Cara berpikir sederhana

Kalau pertanyaannya adalah:

> struktur ini, setelah diberi spectrum gempa, berapa responsnya berdasarkan superposisi modal?

maka itu `SOL111`.

## SOL112

### Makna

Di workflow MYSTRAN kita sekarang, `SOL112` juga masuk ke jalur `MFREQ` dan dipakai sebagai response spectrum run berbasis modal.

Jadi dari sisi ide besar, `SOL112` saat ini bukan dunia yang terpisah total dari `SOL111`. Keduanya sama-sama berada di keluarga:

- modal extraction
- modal participation
- response spectrum per mode
- optional combination

### Kenapa terasa mirip SOL111

Karena secara internal yang sedang kita bangun memang memakai plumbing `MFREQ` yang sama besar.

Jadi untuk user, cara berpikir praktisnya adalah:

- `SOL103` = modal only
- `SOL111` = modal-based response spectrum
- `SOL112` = modal-based response spectrum juga, pada jalur dan workflow deck/output yang sedang kita dukung

### Catatan penting

Kalau nanti arsitektur MYSTRAN dibedakan lebih tegas antara `SOL111` dan `SOL112`, catatan ini bisa direvisi. Tetapi untuk status implementasi sekarang, alurnya memang masih sangat dekat.

## Hasil arah dasar vs combo

Ini bagian yang paling sering bikin bingung.

### Hasil arah dasar

Yang paling dasar adalah hasil per arah:

- `RSX`
- `RSY`
- `RSZ`

Ini sebaiknya dipahami sebagai **hasil primer**.

### Combo

Setelah hasil primer ada, baru boleh dilakukan combination, misalnya:

- `SRSS`
- `CQC`
- orthogonal rule, misalnya `100/30`
- nanti bisa juga percentage rule yang configurable

### Aturan desain yang kita pegang sekarang

- default = **tidak ada combo**
- kalau user tidak minta combo, jangan otomatis bikin hasil gabungan
- cukup tulis hasil arah dasar yang memang diminta

Jadi struktur mentalnya adalah:

```text
modal solve
   -> RSX
   -> RSY
   -> RSZ
   -> optional combination stage
```

bukan:

```text
modal solve
   -> langsung SRSS / CQC tanpa hasil primer yang jelas
```

## Hubungan dengan participation factor dan effective mass

Participation factor dan effective mass itu milik dunia modal, tetapi tetap relevan untuk response spectrum.

Kenapa?

Karena respons spektrum berbasis modal memerlukan:

1. mode shape
2. participation factor tiap mode
3. kontribusi massa efektif tiap mode
4. respons spektrum pada frekuensi/periode mode itu

Jadi tabel-tabel seperti:

- participation factor table
- effective modal mass table

masuk akal muncul pada `SOL103`, dan juga tetap masuk akal muncul pada `SOL111/SOL112`.

Bedanya:

- pada `SOL103`, tabel itu adalah **hasil akhir utama**
- pada `SOL111/SOL112`, tabel itu adalah **hasil antara yang penting**, lalu masih ada tahap respons spektrum di atasnya

## Analogi sederhana

Kalau mau analogi yang mudah:

- `SOL103` = kita hanya memetakan karakter alami struktur
- `SOL111/SOL112` = kita memakai karakter alami itu untuk meramalkan respons terhadap spektrum gempa

Jadi:

- `SOL103` = tahu “siapa dia”
- `SOL111/SOL112` = tahu “bagaimana dia bereaksi”

## Untuk debugging

Kalau hasil response spectrum salah, urutan cek yang sehat biasanya:

1. cek `SOL103` atau modal basis dulu
   - period cocok atau tidak
   - mode shape masuk akal atau tidak
   - participation/effective mass masuk akal atau tidak
2. baru cek input spectrum
3. baru cek directional response `RSX/RSY/RSZ`
4. terakhir cek combo `SRSS/CQC/orthogonal`

Ini penting, karena kalau modal basis sudah salah, hasil RS hampir pasti ikut salah.

## Ringkasan praktis

### Pakai SOL103 ketika

- ingin frekuensi alami
- ingin period
- ingin mode shape
- ingin participation factor
- ingin effective modal mass
- belum ingin response spectrum

### Pakai SOL111 ketika

- ingin response spectrum berbasis modal
- ingin lihat hasil arah dasar lalu optional combo
- ingin tetap punya tabel modal sebagai fondasi perhitungan

### Pakai SOL112 ketika

- ingin response spectrum berbasis modal pada jalur workflow yang sekarang sedang kita hidupkan di MYSTRAN
- ingin output RS yang mengikuti plumbing `MFREQ` aktif kita
- secara konsep besar masih sama keluarga dengan `SOL111`

## Status implementasi lokal saat catatan ini ditulis

Di tree lokal `mystran3`:

- `SOL111` dan `SOL112` keduanya dibawa ke jalur `MFREQ`
- tabel participation factor dan effective modal mass sudah bisa ikut dicetak pada run response spectrum
- hasil primer arah dan combo sekarang sedang dipisahkan dengan aturan:
  - hasil primer dulu
  - combo hanya jika diminta

Itu sebabnya `SOL111` bisa terlihat “modal”, padahal tetap sedang dipakai untuk response spectrum.
