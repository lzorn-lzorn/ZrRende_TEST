# 命名规则

||规则|例|
|----|----|----|
|类名|首字母大写，驼峰命名法|`RenderEngine`|
|类内共有变量|m_变量名, 下划线连接|`m_render_engine`|
|类内私有变量|_p_变量名, 下划线连接|`_p_render_engine`|
|函数名|首字母大写，驼峰命名法|`RenderScene()`|
|普通变量名|小写，下划线连接|`render_engine`|
|常量名|全部大写，下划线连接|`const int MAX_RENDER_ENGINE=1`|
|宏定义|与常量保持一致|`#define MAX_RENDER_ENGINE 100`|
|枚举类型|与常量保持一致|`enum class RenderType`|
|文件名|全部小写，下划线连接|`render_engine.cpp`|

**特例**:如果你命名的实体与已有 C/C++ 实体相似, 可参考现有命名策略.
Eg: `bigopen()`: 函数名, 参照 `open()` 的形式
`bigpos`: `struct` 或 `class`, 参照 pos 的形式
`sparse_hash_map`: STL 型实体; 参照 STL 命名约定
`LONGLONG_MAX`: 常量, 如同 `INT_MAX`

