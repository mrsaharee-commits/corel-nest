# بناء CorelNestEngine.dll عبر GitHub Actions (بدون Visual Studio)
*(English below)*

هذا هو نفس أسلوب المرة السابقة، لكن لمحرّك **CorelNestEngine** (v0.1.1). خبر جيد: هذه النسخة **مكتفية ذاتياً ولا تعتمد على أي مكتبة خارجية** (لا Clipper2)، لذلك يعمل ملف الأوامر مباشرةً ويبني بسرعة.

## ما تحقّقت منه قبل التسليم
جمّعت المحرّك مع اختبار التحقّق `test_harness.cpp` وشغّلته، والنتيجة: **ALL CHECKS PASSED (0 أخطاء)** — لا تداخل بين القطع، وكل قطعة داخل حدود اللوح، والقطعة الأكبر من اللوح تُرفَض بشكل صحيح. وكل الدوال الإحدى عشرة `CNE_*` تُصدَّر بأسماء نظيفة.

## خطوات الاستخدام
1. أنشئ مستودعاً جديداً على github.com (يمكن أن يكون Private).
2. ارفع محتويات هذا المجلد كما هي، **بما في ذلك مجلد `.github`**. (Add file ▸ Upload files ثم اسحب كل الملفات. إن تعذّر رفع المجلد المخفي، أنشئ الملف يدوياً عبر Add file ▸ Create new file بالاسم `.github/workflows/build-dll.yml` والصق محتواه.)
3. افتح تبويب **Actions**؛ سيعمل `build-CorelNestEngine-dll` تلقائياً (أو اضغط **Run workflow**).
4. بعد انتهائه (نحو دقيقة)، افتح التشغيل ونزّل الملف من قسم **Artifacts**: `CorelNestEngine-x64` وبداخله `CorelNestEngine.dll`.

## ملاحظات
- الناتج **64-بت (x64)** لِـ CorelDRAW 2018 فأحدث. لنسخة 32-بت غيّر `-A x64` إلى `-A Win32` في ملف الـ workflow.
- سير العمل يبني الـ DLL ثم **يشغّل الاختبار الذاتي** (ctest) قبل رفع الملف، فإن فشل الاختبار يتوقّف ويُعلمك — ضمان جودة إضافي.
- الـ DLL يُبنى بـ CRT ساكن (‎/MT‎ في `build_win.bat`؛ وفي CMake `MSVC_RUNTIME_LIBRARY=MultiThreaded`) فلا يحتاج ملفات تشغيل إضافية بجانبه — ضعه مباشرة بجوار ملف الـ GMS.
- بديل محلي بدون GitHub: `build_win.bat` (يحتاج "x64 Native Tools Command Prompt" من Visual Studio)، أو استخدم w64devkit مع أمر مماثل (هذه النسخة تُبنى بمترجم واحد لأنها بلا اعتماديات).

---

# Building CorelNestEngine.dll via GitHub Actions (no Visual Studio)

Same approach as before, for the **CorelNestEngine** (v0.1.1) engine. This version is **fully self-contained (no external libraries)**, so the workflow builds quickly with no dependency fetching.

**Verified before delivery:** I built the engine + `test_harness.cpp` and ran it — **ALL CHECKS PASSED (0 failures)**: no overlaps, every part within sheet bounds, oversized part correctly rejected. All 11 `CNE_*` functions export with clean names.

**Steps:** create a GitHub repo, upload this folder (including `.github`), open the **Actions** tab; the `build-CorelNestEngine-dll` workflow builds `CorelNestEngine.dll` on an MSVC Windows runner, runs the self-test, and offers the DLL under the run's **Artifacts** (`CorelNestEngine-x64`).

The workflow also runs the self-test (ctest) before uploading, so a broken build fails loudly. The DLL uses a static CRT, so it needs no extra runtime DLLs beside it. For a 32-bit build, change `-A x64` to `-A Win32`.
