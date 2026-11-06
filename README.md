# Прогноз индекса селективности молекул против VACV

Проект для [Molecular Activity Prediction](https://www.kaggle.com/competitions/molecular-activity-prediction-hack). Задача — предсказать `SI VACV` по SMILES. Официальная метрика — **MSE в исходном масштабе SI**, меньшее значение лучше. Индекс селективности определяется как CC50 / IC50 при сопоставимых измерениях и одинаковых единицах концентрации.

Вход модели ограничен признаками, рассчитанными из структуры молекулы. Колонки `CC50_uM` и `IC50_uM_VACV` запрещены организатором как признаки. Название соединения, штамм, клеточная культура, DOI и источник также не используются для прогноза.

## Состояние проекта

Исходники находятся в `vacv-source.zip`: распакуйте архив в корень проекта перед запуском. Сохранены подготовка данных, загрузка и валидация внешних источников, обучение, прогнозирование и тесты. Полный прогон сравнения ещё выполняется; итоговый submission и результаты пока не утверждены.

Идёт реализация и проверка воспроизводимого решения. Метрика и формат ответа сверены с Kaggle 6 ноября 2026 года. Численные результаты CV, выбранная модель и балл Kaggle добавляются только после фактических запусков. Наличие этого README само по себе не означает, что обучение, проверка в чистом окружении или отправка завершены.

## Данные

Соревнование приватное. Войдите в Kaggle под своим аккаунтом, получите доступ к соревнованию, примите его правила и скачайте файлы во вкладке **Data**. Поместите оригиналы без изменений в `data/raw/`:

```text
data/raw/train.csv
data/raw/test.csv
data/raw/sample_submit.csv
```

В текстовом описании Kaggle образец назван `sample_submission.csv`, однако в списке файлов указан `sample_submit.csv`. Используйте фактически скачанный образец. Формат итогового файла — две колонки в таком порядке: `index`, `SI VACV`. Предсказания должны покрывать все строки test, сохранять исходный порядок идентификаторов и содержать конечные числа.

Исходные файлы соревнования, очищенные копии, фолды со структурами, OOF, прогнозы и обученные веса остаются локальными и исключены через `.gitignore`. Перед публикацией следует отдельно проверить, что текстовые отчёты не включают строки приватного датасета. Исходные материалы задания не копируются в публичный репозиторий.

## Окружение

Начальное окружение: Python 3.9, RDKit 2024.09.6. Точные версии зависимостей зафиксированы в `requirements.txt`.

```bash
python3 -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements.txt
```

Первичная установка требует доступа к реестру Python. Для повторения экспериментов используйте те же версии библиотек и seed 42.

## Порядок работы

Последовательность по ТЗ: загрузка внешних источников → валидация → подготовка датасетов и фолдов → обучение → формирование ответа. Команды выполняются из корня репозитория после установки зависимостей и размещения файлов Kaggle:

```bash
python scripts/download_external_data.py
python scripts/validate_data.py
python scripts/build_dataset.py
python scripts/train.py
python scripts/predict.py
```

Эта последовательность проходит проверку при реализации; завершённость запуска отражается в отчётах. Основной вариант должен воспроизводиться на исходном train независимо от доступности внешних источников. Сводная таблица экспериментов сохраняется в `reports/experiments.csv`, разбор ошибок — в `reports/error_analysis.md`, итоговый ответ — в `submission.csv`. Файл ответа создаётся локально и автоматически не отправляется на Kaggle.

| Каталог | Назначение |
| --- | --- |
| `config/` | Поля входов, метрика, подтверждения правил и конфигурации |
| `vacv/` | Общая стандартизация, признаки и логика моделей |
| `scripts/` | Подготовка данных, обучение и прогнозирование |
| `common/` | Общие фолды и контрольные суммы |
| `external/` | Внешние измерения, кэш и журнал проверки |
| `classical/` | Результаты классических моделей |
| `advanced/` | Результаты дополнительных моделей |
| `reports/` | Аудит, метрики и анализ ошибок |
| `tests/` | Проверки критичных правил и воспроизводимости |
| `docs/` | Разбор задачи и исходных материалов |

## Принципы проверки

Стандартизация по контракту v1: RDKit `MolFromSmiles` с санитизацией → `FragmentParent` → канонический изомерический SMILES и полный InChIKey. Стереохимия сохраняется, дополнительная канонизация таутомеров не применяется. Для каркасов используется Bemis–Murcko без хиральности; все ациклические структуры образуют одну группу `scaffold:ACYCLIC`.

Основная оценка — scaffold CV, контрольная — molecule CV. Группы целиком распределяются между максимум пятью фолдами детерминированно по общему контракту. Одинаковые структуры не должны пересекать обучение и валидацию. Заполнение пропусков, масштабирование и настройка модели выполняются только на обучающей части. Прогноз после обучения на `log1p(SI)` возвращается в исходный масштаб перед оценкой MSE.

Внешние измерения принимаются только после проверки структуры, вируса, условий, единиц и происхождения обоих чисел. Знаки `<` и `>` не превращаются в точные метки. Внешние совпадения с test исключаются; совпадения с валидацией и её каркасами фильтруются отдельно в каждом фолде. Результаты с внешними данными сравниваются с исходным train на одних и тех же фолдах.

## Источники и ограничения

- [Kaggle Overview и Evaluation](https://www.kaggle.com/competitions/molecular-activity-prediction-hack/overview): постановка, метрика и формат ответа.
- [Kaggle Data](https://www.kaggle.com/competitions/molecular-activity-prediction-hack/data): входные файлы, запрет признаков CC50/IC50, требование сбора внешних данных.
- [Kaggle Rules](https://www.kaggle.com/competitions/molecular-activity-prediction-hack/rules): условия участия и обмена кодом/данными.
- Локальные материалы задания перечислены в [обзоре материалов](docs/materials_review.md).

Сбор внешних данных предусмотрен постановкой. Разрешение на внешние предобученные веса отдельно не подтверждено, поэтому они не являются необходимой частью базового решения. Условия каждого внешнего источника фиксируются отдельно; доступность статьи не означает разрешения перепубликовывать весь её текст или данные без проверки лицензии. Публичный репозиторий не заменяет требуемую правилами Kaggle публикацию кода на форуме или в notebook соревнования, если код решено открыть другим участникам во время соревнования.


## Expanded public VACV activity pool

`external/chembl_vacv_activities.csv` contains **4,935** activity records for Vaccinia virus targets CHEMBL613493 and CHEMBL613168, covering **3,594** distinct ChEMBL molecules. Each row includes a canonical SMILES, activity type/value/unit, assay ID, document ID and ChEMBL link. The download is reproducible with `python scripts/download_chembl_vacv.py`; counts by measurement are in `external/chembl_vacv_summary.json`. Source: [ChEMBL activity API](https://www.ebi.ac.uk/chembl/api/data/activity.json?target_chembl_id=CHEMBL613493&limit=1).

**These are not 4,935 extra SI labels.** Only 4 rows are explicitly typed `Selectivity index`; the rest include IC50, EC50, MIC and other endpoints. They are retained as a transparent search pool for matching potency and toxicity under comparable conditions. Do not concatenate the pool with competition training labels. The earlier 12 literature candidates likewise remain under review, not accepted training rows.

## CC50 matched to VACV publications

`external/chembl_cc50_same_molecule_publication.csv` contains **4,457 CC50 measurement rows** for molecules measured against VACV in the **same ChEMBL publication**. The full query found 5,448 CC50 rows across 277 VACV publications; 4,457 match both molecule and document. `external/chembl_cc50_vacv_pair_candidates.csv` lists **340 exact-valued nM potency/CC50 pair candidates** for 192 molecules. See `external/chembl_cc50_summary.json` and reproduce with `python scripts/download_chembl_cc50.py`. [ChEMBL API documentation](https://chembl.gitbook.io/chembl-interface-documentation/web-services/chembl-data-web-services).

Every candidate has `review_status=review_cell_line_assay_and_leakage`. Publication and molecule matching alone do not establish comparable cell lines or assay conditions. **Zero new SI training rows have been accepted** from these pairs. Check those conditions and competition overlap before computing CC50/IC50 (or an appropriate EC50 ratio).

## One public training table

[`vacv_multitask_training.csv`](vacv_multitask_training.csv) is the single public **multi-task training table**: 2,356 rows for 1,201 distinct molecules after excluding structures found in the private competition train/test. It has 1,674 exact CC50, 236 exact VACV IC50, 405 exact VACV EC50 and 41 **provisional** SI rows calculated from exact IC50/CC50 in the same publication and named cell line. Each row carries `endpoint`, `value`, `unit`, assay/document IDs and provenance. See [`vacv_multitask_training_summary.json`](vacv_multitask_training_summary.json).

Train separate output heads or separate models by `endpoint`; do not treat CC50, IC50, EC50 and SI as a single numeric target. SI rows require final assay-method review before being used as definitive labels. For competition SI training, the original private Kaggle train remains local and should be combined only after leakage checks. Rebuild with `python scripts/build_public_multitask_dataset.py --exclude-smiles-csv path/to/train.csv --exclude-smiles-csv path/to/test.csv`. The private files are used only for exclusion and are not published.
