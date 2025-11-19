/*
MIT License

Copyright (c) 2018-2020 Jonathan Young

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
*/
#pragma once
#include <vector>
#include <bx/bx.h>
#include <bx/math.h>
#include <bgfx/bgfx.h>
#include <objzero/objzero.h>
#include <IconFontCppHeaders/IconsFontAwesome4.h>
#include <glm/glm.hpp>
#undef None

namespace ImGui {
	bool Spinner(const char* label);
}

constexpr bgfx::ViewId kModelView = 0;
constexpr bgfx::ViewId kModelTransparentView = 1;
constexpr bgfx::ViewId kGuiView = 2;
constexpr bgfx::ViewId kFirstFreeView = 3;

struct AABB
{
	AABB() : min(FLT_MAX, FLT_MAX, FLT_MAX), max(-FLT_MAX, -FLT_MAX, -FLT_MAX) {}

	void addPoint(bx::Vec3 v)
	{
		min.x = bx::min(min.x, v.x);
		min.y = bx::min(min.y, v.y);
		min.z = bx::min(min.z, v.z);
		max.x = bx::max(max.x, v.x);
		max.y = bx::max(max.y, v.y);
		max.z = bx::max(max.z, v.z);
	}

	void getCorners(bx::Vec3 *corners) const
	{
		const bx::Vec3 aabb[] = { min, max };
		for (int i = 0; i < 8; i++) {
			corners[i].x = aabb[i & 1].x;
			corners[i].y = aabb[(i >> 1) & 1].y;
			corners[i].z = aabb[(i >> 2) & 1].z;
		}
	}

	bx::Vec3 min;
	bx::Vec3 max;
};

void atlasInit();
void atlasShutdown();
void atlasDestroy();
void atlasGenerate();
void atlasFinalize();
void atlasRenderChartsWireframe(const float *modelMatrix);
void atlasShowGuiOptions();
void atlasShowGuiWindow();
uint32_t atlasGetCount();
uint32_t atlasGetWidth();
uint32_t atlasGetHeight();
const uint32_t *atlasGetImage();
struct ModelVertex;
std::vector<ModelVertex> *atlasGetVertices();
std::vector<uint32_t> *atlasGetIndices();
bgfx::VertexBufferHandle atlasGetVb();
bgfx::IndexBufferHandle atlasGetIb();
bgfx::TextureHandle atlasGetFaceDataTexture();
bool atlasIsNotGenerated();
bool atlasIsReady();
void atlasExportChartsImage(const char* filename);
void atlasExportUvData(const char* filename);
void atlasExportRawAtlasImage(const char* filename);
void atlasExportIndividualCharts(const char* baseFilename);
void atlasExportSingleChartImage(uint32_t chartIndex, const char* filename);
void atlasExportSingleChartData(uint32_t chartIndex, const void* chart, const char* filename);

// SPH网格生成器集成
void sphInit();
void sphShutdown();
void sphProcessChart(uint32_t chartIndex);
void sphRunSimulation();
void sphGenerateTriMesh();
void sphGenerateQuadMesh();
void sphGenerateMesh();  // 向后兼容，调用sphGenerateQuadMesh
void sphShowGuiOptions();
// 获取当前选择的边界环数据（用于在Atlas窗口中显示）
const std::vector<glm::vec2>* sphGetCurrentBoundaryLoop();
uint32_t sphGetCurrentBoundaryLoopIndex();
// 获取当前选择的边界环的内环（洞）列表（用于在Atlas窗口中显示）
const std::vector<std::vector<glm::vec2>>* sphGetCurrentBoundaryHoles();
// 获取粒子位置（用于在Atlas窗口中显示）
const std::vector<glm::vec2>* sphGetParticlePositions();
// 批量：处理所有边界环并为每个图表建立模拟
void sphProcessAllBoundaryLoops();
// 批量：运行/单步所有模拟
void sphRunSimulationAll(int steps);
void sphStepSimulationAll();
// 查询：是否已有任意模拟
bool sphHasAnySimulation();
// 获取所有图表的粒子位置（用于3D总览渲染）
const std::vector<glm::vec2>* sphGetAllParticlePositions();
// 控制三维粒子渲染开关
void sphSetShowParticlesIn3D(bool enabled);
bool sphIsShowParticlesIn3D();
// 导出全部图表的归一化边界点（文本）
void sphExportAllBoundaryLoopsNormalized();

// 全局网格显示控制与访问（3D 叠加所有图表的网格）
void sphSetShowAllTriOverlay(bool enabled);
void sphSetShowAllQuadOverlay(bool enabled);
bool sphIsShowAllTriOverlay();
bool sphIsShowAllQuadOverlay();
// 获取所有图表的三/四边形网格（UV 顶点与索引），用于3D渲染
const std::vector<std::vector<glm::vec2>>& sphGetAllTriMeshVertices();
const std::vector<std::vector<glm::uvec3>>& sphGetAllTriMeshIndices();
const std::vector<std::vector<glm::vec2>>& sphGetAllQuadMeshVertices();
const std::vector<std::vector<glm::uvec4>>& sphGetAllQuadMeshIndices();
const std::vector<std::vector<glm::uvec3>>& sphGetAllQuadExtraTriangles();


void sphToggleMeshView();
bool sphIsMeshView();


void sphGetTriMeshForDraw(std::vector<glm::vec2>& outVertices,
                          std::vector<glm::uvec3>& outTriangles);
void sphGetQuadMeshForDraw(std::vector<glm::vec2>& outVertices,
                           std::vector<glm::uvec4>& outQuads);

// Atlas数据访问函数
void* atlasGetData();

void bakeInit();
void bakeShutdown();
void bakeExecute();
void bakeFrame(uint32_t frameNo);
void bakeClear();
void bakeShowGuiOptions();
void bakeShowGuiWindow();
bgfx::TextureHandle bakeGetLightmap();
uint32_t bakeGetLightmapSamplerFlags();
bool bakeIsLightmapReady();
bool bakeIsDenoised();

struct GuiTextureFlags
{
	enum
	{
		PointSampler = 1 << 0,
	};
};

union GuiTexture
{
	void* imgui;
	struct
	{
		bgfx::TextureHandle handle;
		uint16_t flags;
	}
	bgfx;
};

void guiInit();
void guiShutdown();
void guiResize(int width, int height);
void guiRunFrame(float deltaTime);
void guiRender();
bool guiColumnCheckbox(const char *label, const char *id, bool *value);
bool guiColumnColorEdit(const char *label, const char *id, float *color);
bool guiColumnInputFloat(const char *label, const char *id, float *value, float step = 0.0f, float stepFast = 0.0f, const char *format = "%.3f");
bool guiColumnInputInt(const char *label, const char *id, int *value, int step = 1);
bool guiColumnSliderInt(const char *label, const char *id, int *value, int valueMin, int valueMax);

struct ModelVertex
{
	bx::Vec3 pos;
	bx::Vec3 normal;
	float texcoord[4];
	static bgfx::VertexLayout layout;

	static void init()
	{
		layout.begin()
			.add(bgfx::Attrib::Position, 3, bgfx::AttribType::Float)
			.add(bgfx::Attrib::Normal, 3, bgfx::AttribType::Float)
			.add(bgfx::Attrib::TexCoord0, 4, bgfx::AttribType::Float)
			.end();
		assert(layout.getStride() == sizeof(ModelVertex));
	}
};

void modelInit();
void modelShutdown();
void modelFinalize();
void modelOpen(const char *filename);
void modelOpenDialog();
void modelDestroy();
void modelRender(const float *view, const float *projection);
void modelShowGuiMenu();
void modelShowGuiWindow();
AABB modelGetAABB();
const objzModel *modelGetData();
bx::Vec3 modelGetCentroid();
float modelGetScale();
bgfx::ShaderHandle modelGet_vs_model();
bool modelIsLoaded();
bool modelSampleMaterialDiffuse(const objzMaterial *mat, const float *uv, bx::Vec3 *color);
bool modelSampleMaterialEmission(const objzMaterial *mat, const float *uv, bx::Vec3 *color);
// 将2D UV坐标映射到3D位置（用于粒子渲染）
bool modelMapUvTo3D(const float *uv, bx::Vec3 *pos);
// 模型文件名（不含路径与扩展名），用于导出命名
const char* modelGetFileBaseName();

enum class ShadeMode
{
	FlatMaterial,
	LightmapMaterial,
	LightmapOnly
};

enum class WireframeMode
{
	Charts,
	Triangles
};

enum class ChartColorMode
{
	// Sync with xatlas::ChartType
	Planar,
	Ortho,
	LSCM,
	Piecewise,
	Invalid,
	// Not in xatlas::ChartType
	All
};

enum class OverlayMode
{
	None,
	Chart,
	Mesh,
	Stretch
};

struct Options
{
	bool gui = true;
	bool wireframe = true;
	ShadeMode shadeMode = ShadeMode::FlatMaterial;
	WireframeMode wireframeMode = WireframeMode::Triangles;
	ChartColorMode chartColorMode = ChartColorMode::All;
	OverlayMode overlayMode = OverlayMode::None;
	float overlayOpacity = 0.5f;
	int chartCellSize = 1;
	bool lightmapPointSampling = false;
	bool useDenoisedLightmap = true;
	bool showAtlasOptionsWindow = true;
	bool showAtlasWindow = true;
	bool showLightmapWindow = true;
};

extern Options g_options;
struct GLFWwindow;
extern GLFWwindow *g_window;
extern int g_windowSize[2]; // Last known good window size. Not set to 0,0 when minimized.

void randomRGB(uint8_t *color);
uint32_t encodeRGBA(const uint8_t *rgba);
void decodeRGBA(uint32_t rgbaIn, uint8_t *rgbaOut);
void setErrorMessage(const char *format, ...);
void resetCamera();

enum class ShaderId
{
	fs_blit,
	fs_color,
	fs_gui,
	fs_material,
	fs_wireframe,
	vs_blit,
	vs_color,
	vs_gui,
	vs_model,
	vs_wireframe
};

bgfx::ShaderHandle loadShader(ShaderId id);
bgfx::ProgramHandle getColorProgram();
void setWireframeThicknessUniform(float thickness);
bgfx::ProgramHandle getWireframeProgram();

struct WireframeVertex
{
	bx::Vec3 pos;
	bx::Vec3 barycentric;
	static bgfx::VertexLayout layout;

	static void init()
	{
		layout.begin()
			.add(bgfx::Attrib::Position, 3, bgfx::AttribType::Float)
			.add(bgfx::Attrib::TexCoord0, 3, bgfx::AttribType::Float)
			.end();
		assert(layout.getStride() == sizeof(WireframeVertex));
	}
};

#define WINDOW_TITLE "SPHMeshViewer"
