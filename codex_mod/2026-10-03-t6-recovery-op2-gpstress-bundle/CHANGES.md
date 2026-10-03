# Perubahan dari baseline GitHub sampai kondisi akhir

## Baseline dan provenance

Clone awal berada di `C:/PROJECTAI/18a/MYSTRAN`, branch `v18.00.a`, HEAD `931dde4df56d9c39b43d884bda941f7fd27c6036`. Tidak ada commit baru di atas baseline ketika paket ini dibuat. Diff mencakup **seluruh perubahan tracked yang masih ada**, termasuk pekerjaan lokal yang sudah ada sebelum sesi T6 ini.

Folder sibling `C:/PROJECTAI/18a/python` dan `C:/PROJECTAI/18a/test` tidak berada di Git repository tersebut. Oleh sebab itu, snapshot keduanya dimasukkan ke paket sebagai pendukung; hanya file yang memiliki backup awal lokal diberi label perubahan yang dapat dibuktikan. Referensi/dependensi Python yang disalin tidak otomatis berarti file itu telah diedit.

## A. Perubahan lokal yang telah ada ketika pekerjaan ini dilanjutkan

- Kapasitas buffer recovery Q4/CQUADR menggunakan `MELDOF`; assignment BE tertentu dibatasi pada DOF elemen.
- Parameter natural-coordinate MACQ8D dipindah ke scope host pada `CQUAD8_MACQ8D.f90`.
- Output arrays, inisialisasi, dan beberapa pemanggil memakai kapasitas `MAXREQ*5`; perhitungan `MAXREQ` memakai `MAXGROUT + MAXELOUT`.
- Parser/describer menyimpan flag CENTER dan CORNER secara terpisah untuk stress, strain dan force.
- Collector GPSTRESS menerima TRIA6/QUAD8 dan node quadratic.
- Recovery T6 menggunakan displacement basic `UEB`; metadata stress mulai disimpan per elemen.
- `NUM_SEi` T6 berubah dari 4 menjadi 3 dalam data bersama. Output stress/force T6 akhir memakai jumlah eksplisit tujuh titik, sehingga tidak bergantung pada nilai lama ini.

Perubahan tersebut tetap ada dan tercantum dalam patch/daftar file. Tidak ada klaim bahwa semua jalur Q4/Q8, strain atau composite yang terdampak telah tervalidasi numerik.

## B. Referensi Python T6 final

| Selector MYSTRAN | Referensi final |
|---|---|
| SIMOT6 | `Simo1993_Tri6_ShellElement_v2.py` |
| MITC6 | `MITC6_Tri_v4.py` |
| MH6T | `MacNeal_MH6T_Tri_v3.py` |
| REZAIEE | `Rezaiee2017_Tri6_v3.py` |

Registry T6 pada enam benchmark yang ditentukan pengguna diarahkan ke empat referensi di atas. Nama deck `duel3a_mht6_1.dat` diperbaiki menjadi `duel3a_mh6t_1.dat`, dan helper aktif memakai nama yang benar. Registry/formulasi Q8 tidak di-upgrade pada sesi ini.

`compare_f06_vs_python.py` aktif dan `dump_python_debug.py` menggunakan referensi T6 final. Alias lama tetap diterima; label dump memakai ID model seperti E51; urutan natural coordinates mengikuti kelas final. Matriks BM/BB/BS tersedia di CENTER dan keenam node. Koreksi helper juga menghilangkan flip shear global warisan ketika membandingkan kelas final.

## C. Formula T6 diselaraskan dengan Python

- Director geometrik memakai tanda normal centroid yang konsisten. Default geometric directors mengikuti Python; SNORM eksplisit tetap didukung, sedangkan normal hasil averaging otomatis tidak menggantikannya.
- Transformasi covariant-to-local memakai basis pointwise yang sesuai; drilling rotation memakai normal titik.
- MH6T memakai frame centroid geometrik, tanpa pemaksaan sumbu global pada elemen datar. Projection drilling dan basis normal diperbaiki.
- Bending MITC6 dan REZAIEE memakai turunan director serta konvensi cross product final.
- Membrane/shear REZAIEE mengikuti v3: tying locations, rekonstruksi affine, centroid shear dan interpolasi director shear tanpa normalisasi tambahan.

Tujuan tahap ini adalah menyamakan model/solusi dengan Python final sebelum membandingkan recovery. Displacement dibandingkan pada setiap grid dan keenam DOF untuk 2-001, 2-002, 2-003 dan 2-004.

## D. Recovery CENTER, CORNER dan GPSTRESS

Empat kernel kini memulihkan stress/force langsung dalam urutan CENTER lalu node 1–6: `(1/3,1/3)`, `(0,0)`, `(1,0)`, `(0,1)`, `(1/2,0)`, `(1/2,1/2)`, `(0,1/2)`. Quadrature stiffness tetap mengikuti formula masing-masing.

- Semua tujuh sampel disediakan secara internal, termasuk untuk request CENTER, karena format OP2 quadratic memerlukan nilai corner asli.
- F06 mengikuti selector; CORNER T6 juga mencetak tiga node midside.
- Konvensi MYSTRAN menggunakan membrane minus z times curvature dan konversi bending-force negatif. Recovery BE2 menjadi minus Python Bb agar fiber stress dan moment output sesuai; stiffness tidak dibalik tandanya.
- Force recovery menyalin seluruh sampel T6; baris berikut CENTER tidak dibiarkan belum terisi.
- Frame GPSTRESS T6 memakai basis centroid native dengan tanda normal yang konsisten; averaging mencakup enam node, termasuk midsides.
- Index/stride CENTER force T6 mengikuti awal blok elemen yang benar.

STEP A CENTER local, STEP B semua corner local, dan STEP C grid GPSTRESS dibuktikan pada patch 2-001 untuk membrane/bending empat keluarga. Bukti ini tidak otomatis membuktikan recovery pada mesh curved/twisted.

## E. OP2 native T6 dan writer bersama

- Dispatcher OES/OEF mengenali TRIA6; error state T6 dan payload tanpa header OES valid diperbaiki.
- T6 memakai element type 75, bukan label CTRIA3 type 74.
- Stress OES: 70 word/elemen, CENTER + tiga corner, dua fiber; force OEF: 38 word/elemen, CENTER + tiga corner, delapan komponen.
- OES dibuka kembali ketika ITABLE menjadi nol setelah OGS ditutup.
- OGS menulis setiap SURFACE secara terpisah dengan ID asli, menggantikan ID hardcoded 100 yang sebelumnya mencampurkan surface dan menyebabkan ukuran array pyNastran tidak cocok.
- Siemens Nastran 2512 dijalankan pada deck referensi T6. Element type, word count, shape tabel dan ID elemen sesuai dengan MYSTRAN. Ini pemeriksaan kontrak format, bukan klaim semua formula kedua solver identik.
- File campuran `duel3a.OP2` terbaca sampai selesai dengan pyNastran. Nilai OP2 stress/force/GPSTRESS T6 juga dibandingkan dengan F06.

## F. Pembersihan output

Cetakan temporary `RECOV249`, `GPDBG`, norm matriks/director T6 dan `DEBUG_Q8` dihapus. Error diagnostik asli tetap tersedia, begitu juga utilitas Python dump yang diminta pengguna.

Header `FORCES AT GRID POINTS -- SURFACE` kini hanya dicetak jika averaging menghasilkan baris dengan weight positif. Surface kosong untuk keluarga yang sedang ditulis tidak menghasilkan blok kosong. Surface yang memiliki kontribusi nyata tetap tampil. Run ulang `duel3a.dat` memastikan nilai force, GPSTRESS dan displacement tidak berubah.

## Batas validasi dan pekerjaan lanjutan

- Displacement final: 2-001 sampai 2-004, empat keluarga T6, mesh 8/24; regression mesh 2/4 setelah recovery.
- Recovery stress/moment A/B/C: patch 2-001. Stress curved/twisted, thermal, strain/composite dan modal belum ditetapkan oleh gates ini.
- OP2 numerik: T6 OES/OEF dan surface 600 OGS; file mixed structurally terbaca, tetapi kelengkapan numerik setiap Q4/Q8/linear belum diaudit dengan gates ini.
- Direct Femap GUI import belum diuji.
- Kapasitas `MAXREQ*5` yang sudah ada tetap dipertahankan; ini belum dirapikan menjadi alokasi minimal.
- Eksperimen 2-006 dan helper convergence masih memiliki beberapa dependency historis yang hilang. Mereka disertakan sebagai snapshot eksperimen, bukan gate final yang dijamin berjalan. Nama modul tercantum dalam manifest.

## Keputusan satu bundle

Perubahan formula dan output berbagi BE arrays, metadata EID/GID, sizing, frame transforms dan table dispatch. Paket disimpan sebagai satu bundle agar source akhir dan evidence mempunyai kontrak yang konsisten. Pemisahan commit berdasarkan fungsi tetap mungkin sebagai pekerjaan lanjutan, tetapi memerlukan pemilihan hunk dan build/test tiap intermediate commit; hal itu belum dilakukan.
