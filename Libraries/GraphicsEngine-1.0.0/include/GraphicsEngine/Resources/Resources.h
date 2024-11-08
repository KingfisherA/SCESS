#pragma once

namespace GEng {
    class Font;
}

namespace GEng::Res {
    class ShaderSources {
    public:
        static const char *GetColourDiffuseGLSL();
        static const char *GetSolidColourGLSL();
        static const char *GetTextGLSL();
    };

    class FontRes {
    public:
        static Font GetDefaultFont();
    };
}
