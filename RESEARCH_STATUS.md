# الحالة البحثية للمشروع

تاريخ التحديث: **31 أغسطس 2026**.

## الوضع الحالي

هذا المستودع هو حاليًا **مشروع تجريبي/إعادة بناء وقياس**، وليس ادعاءً بامتلاك خوارزمية جديدة لتوليد الأعداد الأولية.

قبل بدء أي اتجاه جديد يجب قراءة:

1. **[PRIOR_ART.md](PRIOR_ART.md)** — ما هو معروف وما الذي قد يبقى مفتوحًا.
2. **[RESULTS_EXTERNAL_BASELINE.md](RESULTS_EXTERNAL_BASELINE.md)** — موقع تطبيقنا أمام `primesieve`.
3. **[RESULTS_BASE_REPRESENTATIONS.md](RESULTS_BASE_REPRESENTATIONS.md)** — دراسة قواعد العد حتى `10^8`.
4. **[RESULTS_COMPACT_CACHE.md](RESULTS_COMPACT_CACHE.md)** — تجارب wheel/cache داخل تطبيقنا.

## Novelty gate

لا نعتبر فكرة جديدة قابلة للبحث قبل الإجابة عن الأسئلة التالية:

- هل توجد في Pritchard / Sorenson / Oliveira e Silva / Járai–Vatai / Atkin–Bernstein / primesieve أو مراجع أخرى؟
- هل هي مجرد إعادة صياغة لـmodular arithmetic؟
- هل يوجد null model واضح إذا كانت النتيجة إحصائية؟
- هل نجحت على validation range مستقل؟
- هل تمت مقارنتها بbaseline خارجي؟
- إذا كان الادعاء أداءً: هل فُصل أثر اللغة والتنفيذ عن أثر الخوارزمية؟

## النتيجة الحالية

- segmented/wheel/primorial/cache-sized sieve: **prior art معروف**.
- تحسينات Python والتمثيل المضغوط: **engineering/replication مفيدة، لا novelty claim**.
- البحث في digit/base statistics بعد إزالة القيود modular: **اتجاه تجريبي مفتوح يحتاج أدبيات أعمق واختبارات إحصائية صارمة**.
- نموذج تنبؤي عام لاختيار wheel/segment/cache: **قد يصبح بحثًا هندسيًا فقط إذا تجاوز heuristics والأعمال السابقة على معماريات متعددة**.

## baseline الحالي

عند `N = 10^8` على نفس GitHub runner:

```text
Project Wheel-30030, Python, enumerate+checksum: ~12.40 s
primesieve C++ iterator, enumerate+checksum:      ~0.0089 s
primesieve single-thread count only:              ~0.006–0.007 s
```

النسبة الكبيرة هنا ليست فرقًا خوارزميًا صرفًا؛ إنها تثبت أن التطبيق الحالي ليس implementation تنافسيًا وأن أي مقارنة خوارزمية جدية تحتاج port منخفض المستوى.
