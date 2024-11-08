#pragma once

#include <cstring>
#include <string>
#include <unordered_map>
#include <glm/matrix.hpp>
#include <vector>

#include "GraphicsEngine/3D/Texture.h"
#include "GraphicsEngine/Utilities/Memory.h"

namespace GEng {
    class VertexBuffer;
    class VertexArray;
    class ShaderProgram;
    class Texture2D;

    struct AtlasGlyphInfo {
        glm::vec2 mPosition;
        glm::vec2 mSize;
        glm::vec2 mBearing;
        float mAdvance;
    };

    using FontAtlasMap = std::unordered_map<char, AtlasGlyphInfo>;

    class Font {
    public:
        explicit Font(const uint8_t *fontData, size_t fontDataSize)
            : mFontDataPtr(fontData), mFontDataSize(fontDataSize) {
        }

        explicit Font(const std::string_view fontPath): mFontPath(fontPath) {
        }

        Font(const Font &other) = default;
        Font(Font &&) noexcept = delete;

        Font &operator=(const Font &other) = default;

        [[nodiscard]] bool IsMemoryFont() const { return mFontDataPtr != nullptr; }

        [[nodiscard]] const uint8_t* GetFontData() const { return mFontDataPtr; }
        [[nodiscard]] size_t GetMemoryFontSize() const { return mFontDataSize; }
        [[nodiscard]] const std::string& GetFontPath() const { return mFontPath; }

    private:
        const uint8_t *mFontDataPtr = nullptr;
        size_t mFontDataSize = 0;
        std::string mFontPath;
    };

    class FontAtlasBase {
    public:
        virtual ~FontAtlasBase() = default;
        [[nodiscard]] virtual const FontAtlasMap &GetAtlasMap() const =0;
        [[nodiscard]] virtual FontAtlasMap &GetAtlasMap() =0;
    };

    class FontAtlasBitmap : public FontAtlasBase {
    public:
        FontAtlasBitmap(): mDimensions(0) {
        }

        FontAtlasBitmap(unsigned char *bitmap, Dimension2D dimensions, const FontAtlasMap &map)
            : mDimensions(dimensions), mAtlasMap(map) {
            size_t noOfBytes = dimensions.x * dimensions.y * sizeof(unsigned char);
            mBitmap = (unsigned char *) malloc(noOfBytes);
            memcpy(mBitmap, bitmap, noOfBytes);
        }

        FontAtlasBitmap(move_construct_t, unsigned char *bitmap, Dimension2D dimensions, FontAtlasMap &&map)
            : mBitmap(bitmap), mDimensions(dimensions), mAtlasMap(std::move(map)) {
        }

        FontAtlasBitmap(const FontAtlasBitmap &other): mDimensions(other.mDimensions), mAtlasMap(other.mAtlasMap) {
            size_t noOfBytes = other.mDimensions.x * other.mDimensions.y * sizeof(unsigned char);
            mBitmap = (unsigned char *) malloc(noOfBytes);
            memcpy(mBitmap, other.mBitmap, noOfBytes);
        }

        FontAtlasBitmap(FontAtlasBitmap &&other) noexcept: mBitmap(other.mBitmap), mDimensions(other.mDimensions),
                                                           mAtlasMap(std::move(other.mAtlasMap)) {
            other.mBitmap = nullptr;
        }

        virtual ~FontAtlasBitmap() {
            free(mBitmap);
        }

        FontAtlasBitmap &operator=(const FontAtlasBitmap &other) {
            mDimensions = other.mDimensions;
            mAtlasMap = other.mAtlasMap;
            size_t noOfBytes = other.mDimensions.x * other.mDimensions.y * sizeof(unsigned char);
            mBitmap = (unsigned char *) malloc(noOfBytes);
            memcpy(mBitmap, other.mBitmap, noOfBytes);

            return *this;
        }

        FontAtlasBitmap &operator=(FontAtlasBitmap &&other) {
            mDimensions = other.mDimensions;
            mAtlasMap = std::move(other.mAtlasMap);
            mBitmap = other.mBitmap;
            other.mBitmap = nullptr;

            return *this;
        }

        virtual const FontAtlasMap &GetAtlasMap() const override { return mAtlasMap; }
        virtual FontAtlasMap &GetAtlasMap() override { return mAtlasMap; }

        const Dimension2D &GetDimensions() const { return mDimensions; }
        Dimension2D &GetDimensions() { return mDimensions; }
        const unsigned char *GetBitmap() const { return mBitmap; }
        unsigned char *&GetBitmap() { return mBitmap; }

    private:
        FontAtlasMap mAtlasMap;
        unsigned char *mBitmap = nullptr;
        Dimension2D mDimensions;
    };

    class FontAtlas : public FontAtlasBase {
    public:
        FontAtlas() = default;

        explicit FontAtlas(const FontAtlasBitmap &bmp): mAtlasMap(bmp.GetAtlasMap()) {
            mTexture = Texture2D::Create(bmp.GetDimensions(), bmp.GetBitmap(), PixelFormat::R);
        }

        FontAtlas(const unsigned char *bitmap, Dimension2D dimensions, const FontAtlasMap &map)
            : mAtlasMap(map) {
            mTexture = Texture2D::Create(dimensions, bitmap, PixelFormat::R);
        }

        FontAtlas(move_construct_t, const unsigned char *bitmap, Dimension2D dimensions, FontAtlasMap &&map)
            : mAtlasMap(std::move(map)) {
            mTexture = Texture2D::Create(dimensions, bitmap, PixelFormat::R);
        }

        virtual const FontAtlasMap &GetAtlasMap() const override { return mAtlasMap; }
        virtual FontAtlasMap &GetAtlasMap() override { return mAtlasMap; }
        const sRef<Texture2D> &GetTexture() const { return mTexture; }

    private:
        FontAtlasMap mAtlasMap;
        sRef<Texture2D> mTexture;
    };

    struct TextLayoutInfo {
    public:
        glm::vec3 bounds;

        enum class Alignment {
            Left = 0,
            Centre,
            Right,
            Justify
        } alignment;
    };

    class FontManager {
    public:
        FontManager();
        ~FontManager();

        enum class GlyphRenderingMode {
            Normal = 0,
            SignedDistanceField
        };


        constexpr static const char *ASCII_GLYPH_SET =
                " !\"#$%&'()*+,-./0123456789:;<=>?@ABCDEFGHIJKLMNOPQRSTUVWXYZ[\\]^_`abcdefghijklmnopqrstuvwxyz{|}~";

        sRef<FontAtlasBitmap> CreateFontAtlasBitmap(const Font& font, const unsigned int glyphResolution,
                                                    const char *glyphSet = ASCII_GLYPH_SET,
                                                    const GlyphRenderingMode mode = GlyphRenderingMode::Normal) const;

        sRef<FontAtlas> CreateFontAtlas(const Font& font, const unsigned int glyphResolution,
                                        const char *glyphSet = ASCII_GLYPH_SET,
                                        const GlyphRenderingMode mode = GlyphRenderingMode::Normal) const;
        static sRef<FontAtlas> CreateFontAtlas(const FontAtlasBitmap &bitmap);

        static const AtlasGlyphInfo &MeasureGlyph(const FontAtlasBase &atlas, const char c) {
            return atlas.GetAtlasMap().at(c);
        }

        // bounds
        // The x and y components define the top left corner of the text run if align is set to left.
        // The z components defines the width of the box. If it is set to zero the text won't be wrapped.
        static std::vector<glm::vec2> LayoutText(const std::string_view text,
                                                 const FontAtlasBase &atlas,
                                                 const TextLayoutInfo &layoutInfo = {
                                                     {0.0f, 0.0f, 0.0f}, TextLayoutInfo::Alignment::Left
                                                 });
        static glm::vec2 MeasureText(const FontAtlasBase &atlas, const std::string_view text);

    private:
        void* mFontLibrary = nullptr;
    };

    class Text {
    public:
        Text();
        ~Text();

        static sRef<Text> Create();

        void SetScaling(const float scaling) { mScaling = scaling; }
        void AddText(const std::string_view text, const TextLayoutInfo &layoutInfo);
        void ClearText();
        void LayoutText(const FontAtlas &atlas);
        void DrawText(const glm::mat4 &pvMatrix) const;

    private:
        std::vector<std::pair<std::string, TextLayoutInfo> > mTextStrings;

        size_t mNoOfGlyphs = 0;
        float mScaling = 0.8f;

        sRef<Texture2D> mFontAtlasTex;
        sRef<VertexBuffer> mTextBuffer;
        sRef<VertexArray> mTextVA;
        static sRef<ShaderProgram> mShader;
    };
}
